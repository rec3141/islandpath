# IslandPath-DIMOB

IslandPath-DIMOB is a standalone software to predict genomic islands in bacterial and archaeal genomes based on the presence of dinucleotide biases and mobility genes.

Genomic islands (GIs) are clusters of genes in prokaryotic genomes of probable horizontal origin.
GIs are disproportionately associated with microbial adaptations of medical or environmental interest.

The latest IslandPath-DIMOB version is integrated in [IslandViewer 4](http://www.pathogenomics.sfu.ca/islandviewer/browse/), the leading integrated web interface for genomic island prediction.

## Python Reimplementation

This fork includes a **Python reimplementation** (`islandpath_dimob.py`) that works directly with GFF3 + FASTA + FAA inputs from tools like Prokka, removing the dependency on GenBank (.gbk) files and broken Perl BioPerl modules.

The Python version implements the same algorithm as the original Perl tool:
1. Dinucleotide relative abundance (4-frame non-overlapping dimer counting); the reference signature pools every CDS of a bin (`--bins`) or, without bins, of a contig
2. 6-gene sliding window bias scoring within each contig (delta* = mean |window_rho - reference_rho|)
3. Thresholding over all windows of the bin or contig (median + 2*SD strong, median + SD moderate; sample SD)
4. Clustering (consecutive marked genes, merge gaps < 6, keep clusters >= 8, drop the last gene)
5. Mobility gene filtering via HMM scan (Pfam gathering thresholds, `--cut_ga`) + keyword matching
6. Islands must contain at least one mobility gene

### Modes and parity with Dimob.pl

`--mode perl-compat` reproduces Dimob.pl (upstream master) island for island:
it gives identical coordinates on *Listeria monocytogenes* EGD-e (NC_003210, 3
islands) and *E. coli* K-12 MG1655 (NC_000913, 17 islands) when GFF3/FASTA/FAA
are generated from the same GenBank record. It exists to test parity.

The default mode departs from Dimob.pl where Dimob.pl is wrong for the biology:

| behaviour | Dimob.pl (perl-compat) | default |
|---|---|---|
| island boundaries | window *i* is labelled with gene *i*+1, so islands shift one gene downstream; marked genes past the last window are dropped, truncating islands at sequence ends | the marked genes, minus the last |
| first gene of a sequence | can never be marked (a 0 sentinel) | can be marked |
| non-ACGT bases | `/[ACGT]{2}/g` restarts after an N, so pairs are double-counted or missed | fixed frame; pairs with non-ACGT skipped |
| HMM mobility hits | full-sequence E-value <= 1e-7: accepts ClpB/ClpC on two weak IstB_IS21 ATPase-domain hits | Pfam gathering thresholds; an IstB hit counts only within 2 genes of another mobility gene (DnaA, DnaC/DnaI also hit IstB) |
| statistics | per genome (Dimob.pl reports only the first record of a multi-record GenBank) | per bin with `--bins`, otherwise per contig |

On MG1655 the default reports 16 islands: 15 are Dimob.pl's islands shifted one
gene upstream; the DnaA-only island at the replication origin is gone; one
island whose only mobility gene is its first gene (b0132, Transposase_31) is
found; and one island whose only mobility evidence is its last gene, FxsA
("suppressor of F exclusion of phage T7", a keyword match), is not.

**Use `--bins` for metagenomes.** Per-contig statistics on a short contig come
from a handful of windows, and a contig that is entirely foreign (a plasmid,
prophage or island fragment) cannot show bias against itself. With
`--bins contig2bin.tsv` (contig<TAB>bin, as DAS Tool writes it) the reference
signature and threshold come from the whole MAG; unbinned contigs fall back to
per-contig statistics.

### Why?

- The original Perl tool requires GenBank (.gbk) input, but Prokka's `tbl2asn` hangs on large metagenomes
- Bioconda's `perl-bioperl 1.7.8` doesn't include `Bio::Perl` (split into `perl-bio-procedural`)
- Metagenome assemblies need per-genome statistics: per bin with `--bins`, otherwise per contig
- The Python version adds a `contig` column to the output (essential for metagenomes)

### Python Dependencies

- Python >= 3.9
- HMMER >= 3.3 (`hmmscan` must be on PATH)
- No other dependencies (pure stdlib)

### Python Usage

```bash
islandpath_dimob.py \
    --gff prokka.gff \
    --fasta assembly.fasta \
    --faa prokka.faa \
    --hmm_db hmmpfam/Pfam-A_mobgenes_201512_prok \
    --cpus 4 \
    --bins contig2bin.tsv \
    -o genomic_islands.tsv
```

Tests: `python3 -m pytest tests`.

Output: tab-separated `island_id`, `contig`, `start`, `end`.

---

## Original Perl Version

### Install

A pre-built Docker image is provided by the [BioContainers project](https://biocontainers.pro/tools/islandpath). Using this pre-installed version of IslandPath-DIMOB ensures the software runs according to expectations.

Users wishing to install locally IslandPath-DIMOB can download a [release](http://github.com/brinkmanlab/islandpath/releases/) and install required perl libraries and HMMER listed below:

Alternatively, you can also clone the latest code from github:

    git clone https://github.com/brinkmanlab/islandpath

Please note that IslandPath-DIMOB predictions should only take a couple of minutes per bacterial genome. It was recently reported that IslandPath-DIMOB was extremely slow on Mac OS X with a conda installation of perl libraries. While we investigate the reason for this issue, we recommend using the Docker image.


**_Dependencies_**

1. Though IslandPath-DIMOB should work with any OS, it has only been tested with linux.

2. Perl version 5.18.x or higher
The latest version of Perl can be obtained from http://www.cpan.org

3. The following Perl libraries are also required:
    - Data::Dumper
    - Log:Log4perl
    - Config::Simple
    - Moose
    - MooseX::Singleton
    - Bio::Perl

4. A working installation of HMMER3
HMMER can be obtained from http://hmmer.org/
"hmmscan" must be within your executable path.
.

### Run

IslandPath-DIMOB v1.0.0 takes as input an annotated complete genome as a genbank (.gbk) or an embl (.embl) file.

    # gbk file
    ./Dimob.pl example/NC_003210.gbk NC_003210_GIs.txt
    # embl file
    ./Dimob.pl example/NC_000913.embl NC_000913_GIs.txt


## Citation

[Bertelli and Brinkman, 2018](https://doi.org/10.1093/bioinformatics/bty095)
[Hsiao et al., 2005](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0010062)


## Questions? Comments? Bugs?

Email islandpick-mail@sfu.ca (contact person: Claire Bertelli) and we'd be happy to help.

If you find a bug, please report it to islandpick-mail@sfu.ca along with the
following information:

* version of Perl (output of 'perl -V' is best)
* version of IslandPath-DIMOB
* operating system type and version
* exact text of error message or description of problem

If we don't have access to a system similar to yours, you may be asked to insert some debugging lines and report back on the results. The more help and information you can provide, the better!


## Copyright and License

IslandPath-DIMOB is distributed under the GNU General Public License. See also the LICENSE file included with this package.


## Versions - New features

### 2026 - Python reimplementation
Python reimplementation of the core algorithm for use with GFF3 + FASTA + FAA inputs.
Removes dependency on GenBank files and broken Perl BioPerl modules.
Adds per-contig processing for metagenome assemblies.

### 23/12/2016 - IslandPath-DIMOB v1.0.0
Increased recall and precision in the prediction of genomic islands based on the presence of dinucleotide bias and mobility genes. Standardization of input file types, and automatic generation of the other file types required by IslandPath-DIMOB.
Input: gbk or embl file
Publication: [Bertelli and Brinkman, 2018](https://doi.org/10.1093/bioinformatics/bty095)

### 2008 - IslandPath-DIMOB
Improvement and assessment of IslandPath-DIMOB predictions by Morgan Langille
Input files: ffn, faa, ptt
Publication: [Langille et al., 2008](http://www.biomedcentral.com/1471-2105/9/329)

### 2005 - IslandPath-DIMOB
Second version developed by Will Hsiao
Further studies used dinucleotide sequence composition bias and the presence of mobility genes to develop a data set of GIs (IslandPath DIMOB) for multiple organisms and revealed that these genomic regions contain higher proportions of novel genes.
Input files: ffn, faa, ptt
Publication: [Hsiao et al., 2005](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0010062)

### 2003 - IslandPath-DINUC
IslandPath-DINUC developed by Will Hsiao
IslandPath was originally designed to aid to the identification of prokaryotic genomics islands (GIs), by visualizing several common characteristics of GIs such as abnormal sequence composition or the presence of genes that functionally related to mobile elements (termed mobility genes).
Publication: [Hsiao et al., 2003](http://bioinformatics.oxfordjournals.org/content/19/3/418.short)


## Authors

IslandPath-DIMOB was written and updated by several members of the [Brinkman Laboratory](http://www.brinkman.mbb.sfu.ca/) at Simon Fraser University, Burnaby, BC, Canada

2015 - present:     Claire Bertelli    claire.bertelli@sfu.ca
2009 - 2015:    Matthew Laird    lairdm@sfu.ca
2007 - 2009: Morgan Langille    morgan.g.i.langille@dal.ca
2003 - 2007: Will Hsiao william.hsiao@bccdc.ca

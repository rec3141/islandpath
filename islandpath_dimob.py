#!/usr/bin/env python3
"""
IslandPath-DIMOB — Genomic island detection via dinucleotide bias + mobility genes.

Python reimplementation of the Perl IslandPath-DIMOB tool (Bertelli et al., 2018;
Hsiao et al., 2003). Works directly with GFF3 + FASTA + FAA outputs from Prokka,
removing the dependency on GenBank (.gbk) input and broken Perl BioPerl modules.

Algorithm:
  1. Group contigs: one group per bin (--bins), otherwise one group per contig
  2. Extract CDS nucleotide sequences from assembly using GFF coordinates
  3. Dinucleotide relative abundance (rho*) via 4-frame non-overlapping counting;
     the reference rho* pools every CDS in the group
  4. 6-gene sliding windows within each contig: delta* = mean |window_rho - group_rho|
  5. Threshold over all windows in the group: strong (median + 2*SD → 6 genes),
     moderate (median + SD → 3 genes); SD is the sample SD, as in Perl
  6. Cluster consecutive marked genes, merge gaps < 6, keep clusters >= 8 genes,
     drop the last gene of each island
  7. Mobility genes: hmmscan with Pfam gathering thresholds (--cut_ga) + keyword
     matching; an IstB (AAA+ ATPase) hit counts only next to another mobility gene
  8. Keep islands with >= 1 mobility gene

--mode perl-compat reproduces Dimob.pl (upstream master) instead, including
behaviour that is wrong for the biology (see PERL_COMPAT below); it exists to
test parity, not for production use.

Reference:
  Bertelli C, et al. (2018) IslandViewer 4. Nucleic Acids Res 46(W1):W59-W64.
  Hsiao W, et al. (2003) IslandPath. Bioinformatics 19(3):418-420.
"""

import argparse
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from urllib.parse import unquote

# ---------------------------------------------------------------------------
# FASTA / GFF parsing
# ---------------------------------------------------------------------------

def parse_fasta(path):
    """Read a FASTA file and return {name: sequence} dict."""
    seqs = {}
    name = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = ''.join(parts)
    return seqs


def parse_gff_cds(path):
    """Extract CDS records per contig from a GFF3 file.

    Returns dict: contig_id -> list of CDS dicts sorted by start position.
    Each CDS dict has keys: contig, start, end, strand, gene_id, product.
    Stops at ##FASTA section (Prokka embeds sequence there).
    """
    cds_by_contig = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            gene_id = ''
            product = ''
            for attr in attrs.split(';'):
                # GFF3 percent-encodes reserved characters in attribute values
                if attr.startswith('ID='):
                    gene_id = unquote(attr[3:])
                elif attr.startswith('product='):
                    product = unquote(attr[8:])
            cds_by_contig[contig].append({
                'contig': contig,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id,
                'product': product,
            })
    # Sort CDS by start position within each contig
    for contig in cds_by_contig:
        cds_by_contig[contig].sort(key=lambda c: c['start'])
    return dict(cds_by_contig)


def revcomp(seq):
    """Reverse complement a DNA sequence."""
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]


def extract_cds_nuc(contig_seq, cds):
    """Extract CDS nucleotide sequence from contig using GFF coordinates.

    GFF is 1-based inclusive. Returns upper-case sequence on the coding strand.
    """
    # GFF coordinates are 1-based, Python slicing is 0-based
    seq = contig_seq[cds['start'] - 1 : cds['end']]
    if cds['strand'] == '-':
        seq = revcomp(seq)
    return seq.upper()

# ---------------------------------------------------------------------------
# Dinucleotide bias calculation
# ---------------------------------------------------------------------------

DINUC_KEYS = [
    'AA', 'AC', 'AG', 'AT',
    'CC', 'CA', 'CG', 'CT',
    'GG', 'GA', 'GC', 'GT',
    'TT', 'TA', 'TC', 'TG',
]

VALID_BASES = set('ACGT')


# Differences between the modes. Each perl-compat entry reproduces Dimob.pl
# master; the default does what the method describes.
#   counting    Perl's SeqWords2 matches /[ACGT]{2}/g, so a non-ACGT base shifts
#               the pair frame: pairs after it are double-counted or missed.
#               Default keeps a fixed frame and skips pairs containing non-ACGT.
#   first gene  Perl seeds its marked list with a 0 sentinel, so gene 0 of a
#               sequence can never be marked. Default can mark it.
#   boundaries  Perl labels window i with gene i+1 and drops marked genes past
#               the last window, shifting islands one gene downstream and
#               truncating islands at sequence ends. Default reports the marked
#               genes, minus the last (the documented precision trade-off).
#   HMM hits    Perl filters on the full-sequence E-value (<= 1e-7), which
#               accepts proteins whose only hits are weak, repeated ATPase
#               domains (e.g. ClpB, ClpC on IstB_IS21). Default uses Pfam's
#               curated gathering thresholds (--cut_ga) and requires an IstB
#               hit to sit next to another mobility gene (IstB acts with IstA;
#               DnaA and DnaI also hit IstB).
#   grouping    Both modes group per contig unless --bins is given; Dimob.pl
#               itself reports only the first record of a multi-record GenBank.
MODES = ('default', 'perl-compat')

_PERL_DIMER_RE = re.compile('[ACGT]{2}')


def count_dimers_nonoverlap(seq, perl_compat=False):
    """Count non-overlapping dimers in seq (step=2).

    Returns (mono_counts, di_counts) as dicts.
    Matches original Perl SeqWords2: step-2 counting, NOT overlapping.
    """
    mono = defaultdict(int)
    di = defaultdict(int)
    # Ensure even length (truncate last base if odd)
    end = len(seq) if len(seq) % 2 == 0 else len(seq) - 1
    if perl_compat:
        seq = seq[:end]
        for b in seq:
            if b in VALID_BASES:
                mono[b] += 1
        for m in _PERL_DIMER_RE.findall(seq):
            di[m] += 1
        return mono, di
    for i in range(0, end, 2):
        b1 = seq[i]
        b2 = seq[i + 1] if i + 1 < len(seq) else None
        if b1 in VALID_BASES:
            mono[b1] += 1
        if b2 is not None and b2 in VALID_BASES:
            mono[b2] += 1
        if b2 is not None and b1 in VALID_BASES and b2 in VALID_BASES:
            di[b1 + b2] += 1
    return mono, di


def count_dimers_4frames(seq, perl_compat=False):
    """Count dimers across 4 frames: original, shifted+1, revcomp, revcomp+shifted.

    This matches the original Perl implementation which combines counts from all
    4 reading frames to get a robust estimate.
    Returns combined (mono_counts, di_counts).
    """
    total_mono = defaultdict(int)
    total_di = defaultdict(int)

    rc = revcomp(seq)

    for s in [seq, seq[1:], rc, rc[1:]]:
        mono, di = count_dimers_nonoverlap(s, perl_compat)
        for k, v in mono.items():
            total_mono[k] += v
        for k, v in di.items():
            total_di[k] += v

    # Ensure all 16 dinucleotide keys are present
    for dk in DINUC_KEYS:
        if dk not in total_di:
            total_di[dk] = 0
    return total_mono, total_di


def calc_relative_abundance(mono, di):
    """Calculate rho*(XY) = f(XY) / (f(X) * f(Y)) for all 16 dinucleotides.

    f(XY) = count(XY) / total_dinucs
    f(X)  = count(X)  / total_mononucs
    """
    total_mono = sum(mono.values())
    total_di = sum(di.values())
    if total_mono == 0 or total_di == 0:
        return {dk: 1.0 for dk in DINUC_KEYS}

    rho = {}
    for dk in DINUC_KEYS:
        b1, b2 = dk[0], dk[1]
        f_di = di.get(dk, 0) / total_di
        f_b1 = mono.get(b1, 0) / total_mono
        f_b2 = mono.get(b2, 0) / total_mono
        denom = f_b1 * f_b2
        rho[dk] = f_di / denom if denom > 0 else 1.0
    return rho


def cds_counts(cds_list, contig_seq, perl_compat=False):
    """Per-CDS (mono, di) counts over 4 frames, in CDS order."""
    counts = []
    for cds in cds_list:
        nuc_seq = extract_cds_nuc(contig_seq, cds)
        counts.append(count_dimers_4frames(nuc_seq, perl_compat))
    return counts


def reference_rho(count_lists):
    """Pooled rho* over every CDS in a group (a genome, bin or contig)."""
    ref_mono = defaultdict(int)
    ref_di = defaultdict(int)
    for counts in count_lists:
        for mono, di in counts:
            for k, v in mono.items():
                ref_mono[k] += v
            for k, v in di.items():
                ref_di[k] += v
    return calc_relative_abundance(ref_mono, ref_di)


def window_biases(counts, ref_rho):
    """delta* (x1000) for each 6-CDS sliding window of one contig.

    Windows never span contigs. Length = max(0, n_cds - 5).
    """
    biases = []
    for i in range(len(counts) - 5):
        win_mono = defaultdict(int)
        win_di = defaultdict(int)
        for mono, di in counts[i:i + 6]:
            for k, v in mono.items():
                win_mono[k] += v
            for k, v in di.items():
                win_di[k] += v
        win_rho = calc_relative_abundance(win_mono, win_di)
        # delta* = mean |window - reference| across 16 dinucleotides, * 1000
        delta = sum(abs(win_rho[dk] - ref_rho[dk]) for dk in DINUC_KEYS) / 16.0
        biases.append(delta * 1000)
    return biases


def bias_threshold_stats(all_biases):
    """(median, sample SD) of the window biases of a group; SD is None if undefined."""
    n = len(all_biases)
    if n < 2:
        return None, None
    v = sorted(all_biases)
    median = v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0
    mean_val = sum(v) / n
    # Statistics::Descriptive::standard_deviation (Perl) is the sample SD
    sd = math.sqrt(sum((x - mean_val) ** 2 for x in v) / (n - 1))
    return median, sd

# ---------------------------------------------------------------------------
# Island detection (thresholding + clustering)
# ---------------------------------------------------------------------------

def find_dinuc_islands(biases, median, sd, min_cluster=8, perl_compat=False):
    """Identify clusters of biased genes from one contig's window biases.

    Thresholds (from original), computed over the whole group:
      - strong: median + 2*SD → mark all 6 genes in the window
      - moderate: median + SD → mark first 3 genes in the window

    Clustering:
      - Collect marked gene indices
      - Find consecutive runs
      - Keep runs >= min_cluster genes
      - Merge runs with gap < 6 genes apart
      - Exclude last gene of each island (per original precision tuning)

    Returns list of (start_idx, end_idx) tuples into the CDS list (inclusive).
    """
    if not biases or not sd:
        return []

    strong_thresh = median + 2 * sd
    moderate_thresh = median + sd

    # Mark gene indices that show bias
    marked = set()
    for i, bias in enumerate(biases):
        if bias > strong_thresh:
            # Mark all 6 genes in this window
            for j in range(6):
                marked.add(i + j)
        elif bias > moderate_thresh:
            # Mark first 3 genes
            for j in range(3):
                marked.add(i + j)
    if perl_compat:
        marked.discard(0)

    if not marked:
        return []

    # Sort and find consecutive clusters
    marked_sorted = sorted(marked)

    clusters = []
    current = [marked_sorted[0]]
    for idx in marked_sorted[1:]:
        if idx == current[-1] + 1:
            current.append(idx)
        else:
            clusters.append(current)
            current = [idx]
    clusters.append(current)

    # Keep clusters >= min_cluster genes
    clusters = [c for c in clusters if len(c) >= min_cluster]

    if not clusters:
        return []

    # Merge clusters that are < 6 genes apart
    merged = [clusters[0]]
    for cluster in clusters[1:]:
        gap = cluster[0] - merged[-1][-1] - 1
        if gap < 6:
            # Fill the gap
            bridge = list(range(merged[-1][-1] + 1, cluster[0]))
            merged[-1] = merged[-1] + bridge + cluster
        else:
            merged.append(cluster)

    # Exclude last gene of each island (per original: improves precision)
    islands = []
    for cluster in merged:
        if perl_compat:
            # Perl resolves marked index k through window k's label, which
            # names gene k+1; indices past the last window have no label.
            kept = [k for k in cluster[:-1] if k < len(biases)]
            if kept:
                islands.append((kept[0] + 1, kept[-1] + 1))
        elif len(cluster) > 1:
            islands.append((cluster[0], cluster[-2]))

    return islands

# ---------------------------------------------------------------------------
# Mobility gene detection
# ---------------------------------------------------------------------------

MOBILITY_KEYWORDS = [
    'transposase',
    'istb',
    'insertion element',
    'recombinase',
    'insertion sequence',
    'resolvase',
    'integrase',
    'phage',
    'transposon',
    'transposable element',
    'excisionase',
]

MOBILITY_PATTERNS = [re.compile(kw, re.IGNORECASE) for kw in MOBILITY_KEYWORDS]


def find_keyword_mobility(cds_list):
    """Return set of gene_ids whose product annotation matches mobility keywords."""
    mob_ids = set()
    for cds in cds_list:
        product = cds.get('product', '')
        if any(pat.search(product) for pat in MOBILITY_PATTERNS):
            mob_ids.add(cds['gene_id'])
    return mob_ids


ISTB_MODELS = {'IstB_IS21', 'IstB_IS21_ATP'}


def run_hmmscan(faa_path, hmm_db, cpus=1, perl_compat=False, evalue=1e-7):
    """Run hmmscan; return {protein ID: set of mobility HMM names hit}.

    Default: Pfam gathering thresholds (--cut_ga) decide what counts as a hit.
    perl-compat: full-sequence E-value <= evalue (Bio::SearchIO significance).
    """
    hits = defaultdict(set)

    if not os.path.isfile(faa_path) or os.path.getsize(faa_path) == 0:
        return hits

    with tempfile.NamedTemporaryFile(mode='w', suffix='.domtbl', delete=False) as tmp:
        domtbl_path = tmp.name

    try:
        cmd = ['hmmscan', '--domtblout', domtbl_path, '--noali', '--cpu', str(cpus)]
        cmd += ['-E', str(evalue)] if perl_compat else ['--cut_ga']
        cmd += [hmm_db, faa_path]
        res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            # A failed scan must not look like "no mobility genes".
            raise RuntimeError(f"hmmscan failed (exit {res.returncode}): {res.stderr.strip()}")
        with open(domtbl_path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) < 13:
                    continue
                model, query_name = fields[0], fields[3]
                # perl-compat: full-sequence E-value (col 7); with --cut_ga,
                # hmmscan reports only hits that pass the curated thresholds.
                if not perl_compat or float(fields[6]) <= evalue:
                    hits[query_name].add(model)
    finally:
        if os.path.isfile(domtbl_path):
            os.unlink(domtbl_path)

    return hits


def mobility_gene_ids(cds_list, hmm_hits, perl_compat=False, istb_neighbourhood=2):
    """IDs of mobility genes on one contig (HMM hits + product keywords).

    Default: a gene whose only evidence is an IstB HMM hit counts only if
    another mobility gene lies within istb_neighbourhood genes on the contig.
    """
    keyword = find_keyword_mobility(cds_list)
    ids = [c['gene_id'] for c in cds_list]
    if perl_compat:
        return keyword | {g for g in ids if hmm_hits.get(g)}
    strong = [g in keyword or bool(hmm_hits.get(g, set()) - ISTB_MODELS) for g in ids]
    mob = set()
    for i, g in enumerate(ids):
        if strong[i]:
            mob.add(g)
        elif hmm_hits.get(g):
            lo, hi = max(0, i - istb_neighbourhood), i + istb_neighbourhood + 1
            if any(strong[lo:i] + strong[i + 1:hi]):
                mob.add(g)
    return mob

# ---------------------------------------------------------------------------
# Per-contig pipeline
# ---------------------------------------------------------------------------

def process_group(contigs, assembly, cds_by_contig, hmm_hits, perl_compat=False,
                  min_island_bp=0, min_cds=8):
    """Run IslandPath-DIMOB on one group of contigs (a bin, or a single contig).

    The reference signature and the bias threshold come from the whole group;
    windows and islands stay within contigs. Returns [{contig, start, end}].
    """
    members = [c for c in contigs if len(cds_by_contig.get(c, [])) >= 6]
    if sum(len(cds_by_contig[c]) for c in members) < min_cds:
        return []
    counts = {c: cds_counts(cds_by_contig[c], assembly[c], perl_compat) for c in members}
    ref = reference_rho(counts.values())
    biases = {c: window_biases(counts[c], ref) for c in members}
    median, sd = bias_threshold_stats([b for c in members for b in biases[c]])
    if sd is None:
        return []

    islands = []
    for c in members:
        cds_list = cds_by_contig[c]
        if len(cds_list) < 8:
            continue
        dinuc_islands = find_dinuc_islands(biases[c], median, sd, 8, perl_compat)
        if not dinuc_islands:
            continue
        mob_ids = mobility_gene_ids(cds_list, hmm_hits, perl_compat)
        for start_idx, end_idx in dinuc_islands:
            island_cds = cds_list[start_idx:end_idx + 1]
            if not {x['gene_id'] for x in island_cds} & mob_ids:
                continue
            island_start = island_cds[0]['start']
            island_end = island_cds[-1]['end']
            if island_end - island_start + 1 < min_island_bp:
                continue
            islands.append({'contig': c, 'start': island_start, 'end': island_end})
    return islands


def read_bins(path, assembly):
    """contig -> bin from a two-column TSV (contig<TAB>bin, as DAS Tool writes)."""
    bins, unknown = {}, 0
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            contig, b = line.rstrip('\n').split('\t')[:2]
            if contig in assembly:
                bins[contig] = b
            else:
                unknown += 1
    if not bins:
        raise SystemExit(f'[ERROR] no contig in {path} matches the assembly '
                         '(expected contig<TAB>bin)')
    if unknown:
        print(f'[WARNING] {unknown} lines in {path} name contigs not in the assembly',
              file=sys.stderr)
    return bins

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='IslandPath-DIMOB: genomic island detection via dinucleotide bias + mobility genes',
    )
    parser.add_argument('--gff', required=True, help='Prokka GFF3 file')
    parser.add_argument('--fasta', required=True, help='Assembly FASTA')
    parser.add_argument('--faa', required=True, help='Prokka protein FASTA (.faa)')
    parser.add_argument('--hmm_db', required=True, help='Path to HMM database (Pfam mobility gene profiles)')
    parser.add_argument('--cpus', type=int, default=1, help='CPUs for hmmscan')
    parser.add_argument('--bins', help='contig<TAB>bin table: reference signature and '
                        'threshold per bin instead of per contig (recommended for MAGs)')
    parser.add_argument('--mode', choices=MODES, default='default',
                        help='perl-compat reproduces Dimob.pl, including its bugs '
                        '(for parity testing)')
    parser.add_argument('--min_island_bp', type=int, default=0,
                        help='drop islands shorter than this (Dimob.pl applies none)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    # Validate inputs
    for path, label in [(args.gff, 'GFF'), (args.fasta, 'FASTA'), (args.faa, 'FAA')]:
        if not os.path.isfile(path):
            print(f'[ERROR] {label} file not found: {path}', file=sys.stderr)
            sys.exit(1)

    # Parse inputs
    print(f'[INFO] Parsing assembly FASTA: {args.fasta}', file=sys.stderr)
    assembly = parse_fasta(args.fasta)

    print(f'[INFO] Parsing GFF3: {args.gff}', file=sys.stderr)
    cds_by_contig = parse_gff_cds(args.gff)

    # Run batch HMM scan on all proteins
    print(f'[INFO] Running hmmscan on {args.faa} ({args.cpus} CPUs)', file=sys.stderr)
    perl_compat = args.mode == 'perl-compat'
    hmm_hits = run_hmmscan(args.faa, args.hmm_db, cpus=args.cpus, perl_compat=perl_compat)
    print(f'[INFO] Proteins with mobility HMM hits: {len(hmm_hits)}', file=sys.stderr)

    contigs = []
    for contig_id in sorted(cds_by_contig.keys()):
        if contig_id not in assembly:
            print(f'[WARNING] Contig {contig_id} in GFF but not in FASTA — skipping', file=sys.stderr)
            continue
        contigs.append(contig_id)

    # One group per bin; unbinned contigs (or every contig, without --bins) alone
    groups = defaultdict(list)
    bins = read_bins(args.bins, assembly) if args.bins else {}
    for c in contigs:
        groups[('bin', bins[c]) if c in bins else ('contig', c)].append(c)
    print(f'[INFO] Mode: {args.mode}; groups: {sum(k[0] == "bin" for k in groups)} bins, '
          f'{sum(k[0] == "contig" for k in groups)} single contigs', file=sys.stderr)

    all_islands = []
    for key in sorted(groups):
        all_islands.extend(process_group(groups[key], assembly, cds_by_contig, hmm_hits,
                                         perl_compat, args.min_island_bp))
    all_islands.sort(key=lambda x: (x['contig'], x['start']))
    contigs_processed = len(contigs)
    contigs_with_islands = len({x['contig'] for x in all_islands})

    # Write output
    print(f'[INFO] Contigs processed: {contigs_processed}', file=sys.stderr)
    print(f'[INFO] Contigs with islands: {contigs_with_islands}', file=sys.stderr)
    print(f'[INFO] Total islands found: {len(all_islands)}', file=sys.stderr)

    with open(args.output, 'w') as fh:
        fh.write('island_id\tcontig\tstart\tend\n')
        for i, island in enumerate(all_islands, 1):
            fh.write(f"GI_{i}\t{island['contig']}\t{island['start']}\t{island['end']}\n")

    print(f'[INFO] Results written to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()

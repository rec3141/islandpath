"""Unit tests for the default / perl-compat differences (run: python3 -m pytest tests)."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import islandpath_dimob as ip  # noqa: E402


def test_non_acgt_base_does_not_shift_the_pair_frame():
    # The shifted frame of ACGNTTA is CGNTTA. A fixed frame gives CG, (NT
    # skipped), TA. Perl's /[ACGT]{2}/g restarts after the N and gives CG, TT,
    # so TT is counted in both frames and TA in neither.
    _, di = ip.count_dimers_nonoverlap('ACGNTT')
    assert dict(di) == {'AC': 1, 'TT': 1}
    _, di = ip.count_dimers_nonoverlap('CGNTTA')
    assert dict(di) == {'CG': 1, 'TA': 1}
    _, di = ip.count_dimers_nonoverlap('CGNTTA', perl_compat=True)
    assert dict(di) == {'CG': 1, 'TT': 1}


def _strong_run(n_windows, start, length):
    biases = [0.0] * n_windows
    for i in range(start, start + length):
        biases[i] = 100.0
    return biases


def test_first_gene_can_be_marked_only_in_default():
    biases = _strong_run(40, 0, 4)          # windows 0-3 -> genes 0-8
    median, sd = 0.0, 1.0
    assert ip.find_dinuc_islands(biases, median, sd) == [(0, 7)]
    # Perl never marks gene 0 and shifts by one gene: genes 2..8
    assert ip.find_dinuc_islands(biases, median, sd, perl_compat=True) == [(2, 8)]


def test_perl_compat_truncates_islands_at_sequence_end():
    biases = _strong_run(20, 14, 6)         # windows 14-19 -> genes 14-24 (last gene)
    median, sd = 0.0, 1.0
    assert ip.find_dinuc_islands(biases, median, sd) == [(14, 23)]
    # Marked genes past the last window (index 19) have no label in Perl
    assert ip.find_dinuc_islands(biases, median, sd, perl_compat=True) == [(15, 20)]


def test_sample_sd():
    median, sd = ip.bias_threshold_stats([1.0, 2.0, 3.0, 4.0])
    assert median == 2.5
    assert sd == pytest.approx(1.2909944, rel=1e-6)


def _cds(names, products=None):
    products = products or {}
    return [{'gene_id': g, 'product': products.get(g, '')} for g in names]


def test_istb_needs_a_mobility_neighbour_in_default():
    cds = _cds(['dnaA', 'g1', 'g2', 'g3', 'istB', 'istA'])
    hits = {'dnaA': {'IstB_IS21'}, 'istB': {'IstB_IS21'}, 'istA': {'rve'}}
    assert ip.mobility_gene_ids(cds, hits) == {'istB', 'istA'}
    assert ip.mobility_gene_ids(cds, hits, perl_compat=True) == {'dnaA', 'istB', 'istA'}


def test_keyword_neighbour_rescues_istb():
    cds = _cds(['istB', 'tnp'], {'tnp': 'IS21 family transposase'})
    assert ip.mobility_gene_ids(cds, {'istB': {'IstB_IS21_ATP'}}) == {'istB', 'tnp'}


def test_failed_hmmscan_raises(tmp_path):
    faa = tmp_path / 'p.faa'
    faa.write_text('>p1\nMKV\n')
    with pytest.raises((RuntimeError, FileNotFoundError)):
        ip.run_hmmscan(str(faa), str(tmp_path / 'missing.hmm'))


def test_bins_with_swapped_columns_fail(tmp_path):
    t = tmp_path / 'bins.tsv'
    t.write_text('bin1\tcontig_1\n')
    with pytest.raises(SystemExit):
        ip.read_bins(str(t), {'contig_1': 'ACGT'})


def test_gff_attributes_are_percent_decoded(tmp_path):
    g = tmp_path / 'a.gff'
    g.write_text('c1\tx\tCDS\t1\t9\t.\t+\t0\tID=c1_1;product=lysine decarboxylase%2C acid-inducible\n')
    assert ip.parse_gff_cds(str(g))['c1'][0]['product'] == 'lysine decarboxylase, acid-inducible'

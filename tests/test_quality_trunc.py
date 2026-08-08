"""Truncation length selection must not be fooled by fragment reads.

The failure this guards is silent: a bimodal read-length distribution drags the
10th percentile into the adapter-dimer mode, truncation lands far below the
amplicon, ~99% of pairs fail to merge, and the sample simply disappears from the
final table with no error anywhere.
"""
import importlib.util
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "microscape_quality", os.path.join(_HERE, "..", "microscape", "quality.py"))
quality = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(quality)


def _bimodal(n_fragment=325, n_real=1101, fragment_len=31, real_len=232):
    """A real 18S library's shape: 22.8% fragments against a 232 bp mode."""
    return [fragment_len] * n_fragment + [real_len] * n_real


def test_fragments_are_excluded_from_the_length_distribution():
    lengths = _bimodal()
    usable = quality._amplicon_lengths(lengths)
    assert len(usable) == 1101
    assert min(usable) == 232


def test_percentile_lands_on_the_amplicon_not_the_fragments():
    # The whole point: p10 over the raw distribution is in the fragment mode.
    lengths = _bimodal()
    assert int(np.percentile(lengths, 10)) < 50
    assert int(np.percentile(quality._amplicon_lengths(lengths), 10)) > 200


def test_clean_distributions_are_left_alone():
    # A healthy library must not be altered — this fix has to be a no-op there.
    lengths = [233] * 30000 + [248] * 481
    assert quality._amplicon_lengths(lengths) == lengths


def test_a_library_that_is_almost_all_fragments_is_not_over_filtered():
    # Not bimodal, just broken. Measuring the 2% remainder would be less honest
    # than measuring everything, so keep the raw distribution.
    lengths = [30] * 990 + [232] * 10
    assert quality._amplicon_lengths(lengths) == lengths


def test_empty_input_is_handled():
    assert quality._amplicon_lengths([]) == []


def test_trunc_pos_recovers_a_usable_length_on_a_bimodal_sample():
    """End to end through _find_trunc_pos, with quality held constant.

    Quality is uniformly high so the length cap is what decides, which is the
    path that used to collapse to the floor.
    """
    quals = [[35] * n for n in _bimodal()]
    pos = quality._find_trunc_pos(quals, min_q=25.0, window=10, min_length=150)
    # Must clear the old floor of 150 by a wide margin: a 277 bp amplicon needs
    # well over 138 bp per read to leave any usable overlap.
    assert pos > 200, f"truncation collapsed to {pos}"


def test_the_floor_never_asks_for_reads_longer_than_the_library():
    """A 2x150 run whose amplicon reads are 126bp after primer removal.

    dada2 discards reads shorter than truncLen, so a 150 floor applied to the
    length cap truncates at 148 (the longest read present) and keeps 607 of
    101194 reads — the sample survives the filter as a rounding error and the
    run reports "no results". The floor belongs on the quality position, not on
    the cap: it exists to stop a quality dip truncating below overlap, not to
    ask for bases that were never sequenced.
    """
    lengths = [126] * 100215 + [124] * 233 + [125] * 33 + [127] * 82 + [148] * 607
    quals = [[35] * n for n in lengths]
    pos = quality._find_trunc_pos(quals, min_q=25.0, window=10, min_length=150)
    assert pos <= 126, f"truncation {pos} exceeds the reads that exist"
    kept = sum(1 for n in lengths if n >= pos)
    assert kept > 0.99 * len(lengths), f"only {kept} of {len(lengths)} reads survive"


def test_the_floor_still_holds_against_a_quality_cliff():
    """Issue #4 must survive the change: a dip at ~30bp still floors to 150."""
    lengths = [300] * 5000
    quals = [[35] * 30 + [5] * 270 for _ in lengths]
    pos = quality._find_trunc_pos(quals, min_q=25.0, window=10, min_length=150)
    assert pos == 150, f"quality cliff truncated to {pos}, not the 150 floor"

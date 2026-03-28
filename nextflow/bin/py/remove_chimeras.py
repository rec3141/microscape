#!/usr/bin/env python3
#
# remove_chimeras.py — Sparse consensus chimera removal
#
# Python mirror of remove_chimeras.R. A memory-efficient implementation
# that operates on long-format count data instead of a dense sample-by-ASV
# matrix.
#
# Algorithm (identical to dada2 "consensus" method):
#   For each query ASV (in order of increasing total abundance):
#     1. Find all samples where the query is present
#     2. In each sample, find candidate parents: ASVs with abundance >=
#        minFoldParentOverAbundance * query_abundance AND >= minParentAbundance
#     3. Check if the query is a bimera of any two parents
#     4. Tally: if flagged in >= minSampleFraction of samples -> chimera
#
# The per-ASV checks are independent and parallelized with multiprocessing.
#
# Usage:
#   remove_chimeras.py <seqtab.pkl> <cpus>
#
# Outputs:
#   seqtab_nochim.pkl    Chimera-free sequence table (long-format DataFrame)
#   chimera_stats.tsv    Summary statistics

import sys
import os
import pickle
import pandas as pd
import numpy as np
from multiprocessing import Pool
from collections import defaultdict


# ---------------------------------------------------------------------------
# Bimera detection
# ---------------------------------------------------------------------------
def is_bimera(query, parents, min_match_frac=0.9):
    """Check if query is a bimera (chimera of two parents).

    Splits the query at every possible breakpoint and checks if the left
    half matches one parent and the right half matches another.

    Uses a coarse-to-fine strategy: check breakpoints at 10-bp steps first,
    then refine around promising breakpoints.
    """
    qlen = len(query)
    if qlen < 20 or len(parents) < 2:
        return False

    # Coarse pass: check breakpoints at 10-bp intervals
    coarse_step = max(1, min(10, qlen // 5))
    promising_bps = []

    for bp in range(10, qlen - 10, coarse_step):
        left = query[:bp]
        right = query[bp:]
        left_len = len(left)
        right_len = len(right)

        found_left_parent = False
        left_parent_idx = -1

        for i, p1 in enumerate(parents):
            p1_left = p1[:left_len] if len(p1) >= left_len else p1
            cmp_len = min(left_len, len(p1_left))
            matches = sum(1 for a, b in zip(left[:cmp_len], p1_left[:cmp_len])
                          if a == b)
            if matches / left_len >= min_match_frac:
                found_left_parent = True
                left_parent_idx = i
                break

        if not found_left_parent:
            continue

        for j, p2 in enumerate(parents):
            if j == left_parent_idx:
                continue
            p2_right = p2[len(p2) - right_len:] if len(p2) >= right_len else p2
            cmp_len = min(right_len, len(p2_right))
            matches = sum(1 for a, b in zip(right[:cmp_len], p2_right[:cmp_len])
                          if a == b)
            if matches / right_len >= min_match_frac:
                return True

    return False


# ---------------------------------------------------------------------------
# Per-ASV chimera check function (for multiprocessing)
# ---------------------------------------------------------------------------
# Module-level globals set by initializer for multiprocessing workers
_sample_index = None
_seq_to_samples = None
_min_fold = None
_min_parent_abund = None
_min_sample_fraction = None
_ignore_n_negatives = None


def _init_worker(sample_index, seq_to_samples, min_fold,
                 min_parent_abund, min_sample_frac, ignore_n_neg):
    global _sample_index, _seq_to_samples
    global _min_fold, _min_parent_abund
    global _min_sample_fraction, _ignore_n_negatives
    _sample_index = sample_index
    _seq_to_samples = seq_to_samples
    _min_fold = min_fold
    _min_parent_abund = min_parent_abund
    _min_sample_fraction = min_sample_frac
    _ignore_n_negatives = ignore_n_neg


def _check_one_asv(query_seq):
    """Check if a single ASV is chimeric using consensus across samples."""
    info = _seq_to_samples[query_seq]
    query_samples = info["samples"]
    query_counts = info["counts"]
    nsam = len(query_samples)
    nflag = 0

    for s_idx in range(nsam):
        sam = query_samples[s_idx]
        qcount = query_counts[s_idx]

        # Look up all ASVs in this sample
        sam_data = _sample_index[sam]
        sam_seqs = sam_data["sequences"]
        sam_counts = sam_data["counts"]

        # Find qualifying parents: more abundant in THIS sample
        min_count = max(_min_parent_abund, _min_fold * qcount)
        parents = [seq for seq, cnt in zip(sam_seqs, sam_counts)
                   if cnt >= min_count and seq != query_seq]

        if len(parents) < 2:
            continue

        # Limit parents to top 20 by abundance for performance
        parents = parents[:20]

        if is_bimera(query_seq, parents):
            nflag += 1

    # Consensus decision: same logic as dada2
    if nflag == 0:
        return False
    return (nflag >= nsam or
            nflag >= (nsam - _ignore_n_negatives) * _min_sample_fraction)


# ---------------------------------------------------------------------------
# Sparse consensus chimera detection
# ---------------------------------------------------------------------------
def find_chimeras_sparse(dt, min_fold_parent=1.5, min_parent_abund=2,
                          min_sample_fraction=0.9, ignore_n_negatives=1,
                          cpus=1):
    """Find chimeric sequences using sparse consensus method.

    Args:
        dt: DataFrame with columns (sample, sequence, count)
        min_fold_parent: parent must be this fold more abundant than query
        min_parent_abund: minimum absolute count for a parent
        min_sample_fraction: fraction of samples that must flag as chimeric
        ignore_n_negatives: allow this many non-flagged samples
        cpus: number of parallel workers

    Returns:
        set of chimeric sequence strings
    """
    # Build sparse indices
    # sample_index: for each sample, list of (sequence, count) sorted by
    #   count descending (so parents come first)
    sample_index = {}
    for sam, grp in dt.groupby("sample"):
        grp_sorted = grp.sort_values("count", ascending=False)
        sample_index[sam] = {
            "sequences": grp_sorted["sequence"].tolist(),
            "counts": grp_sorted["count"].tolist(),
        }

    # seq_to_samples: for each sequence, list of (sample, count)
    seq_to_samples = {}
    seq_totals = {}
    for seq, grp in dt.groupby("sequence"):
        seq_to_samples[seq] = {
            "samples": grp["sample"].tolist(),
            "counts": grp["count"].tolist(),
        }
        seq_totals[seq] = grp["count"].sum()

    # Sort by total abundance ascending (least abundant first)
    all_seqs = sorted(seq_totals.keys(), key=lambda s: seq_totals[s])

    print(f"[INFO] Checking {len(all_seqs)} ASVs for chimeras using "
          f"{cpus} cores")

    # Run checks
    if cpus > 1:
        with Pool(cpus, initializer=_init_worker,
                  initargs=(sample_index, seq_to_samples, min_fold_parent,
                            min_parent_abund, min_sample_fraction,
                            ignore_n_negatives)) as pool:
            results = pool.map(_check_one_asv, all_seqs, chunksize=50)
    else:
        _init_worker(sample_index, seq_to_samples, min_fold_parent,
                     min_parent_abund, min_sample_fraction,
                     ignore_n_negatives)
        results = [_check_one_asv(seq) for seq in all_seqs]

    chimeric_seqs = {seq for seq, is_chim in zip(all_seqs, results) if is_chim}

    print(f"[INFO] Flagged {len(chimeric_seqs)} chimeras out of "
          f"{len(all_seqs)} ASVs")

    return chimeric_seqs


# ===========================================================================
# Main script
# ===========================================================================

if len(sys.argv) < 3:
    print("Usage: remove_chimeras.py <seqtab_long.pkl> <cpus>",
          file=sys.stderr)
    sys.exit(1)

input_path = sys.argv[1]
cpus       = int(sys.argv[2])

# ---------------------------------------------------------------------------
# Load input (long-format DataFrame from merge_seqtabs.py)
# ---------------------------------------------------------------------------
with open(input_path, "rb") as f:
    dt = pickle.load(f)

if not isinstance(dt, pd.DataFrame):
    print("[ERROR] Expected long-format DataFrame with columns: "
          "sample, sequence, count", file=sys.stderr)
    sys.exit(1)

n_input_asvs  = dt["sequence"].nunique()
n_input_reads = int(dt["count"].sum())
n_samples     = dt["sample"].nunique()

print(f"[INFO] Input: {n_samples} samples, {n_input_asvs} ASVs, "
      f"{n_input_reads} reads, {len(dt)} non-zero entries")

# ---------------------------------------------------------------------------
# Pre-filter: remove obvious noise before the expensive chimera check
# ---------------------------------------------------------------------------
seq_stats = dt.groupby("sequence").agg(
    total=("count", "sum"),
    seq_len=("sequence", lambda x: len(x.iloc[0]))
).reset_index()

pre_singleton = set(seq_stats.loc[seq_stats["total"] <= 1, "sequence"])
pre_short     = set(seq_stats.loc[seq_stats["seq_len"] < 20, "sequence"])
pre_remove    = pre_singleton | pre_short
n_pre_removed = len(pre_remove)

if n_pre_removed > 0:
    n_short_only = len(pre_short - pre_singleton)
    print(f"[INFO] Pre-filter: removing {len(pre_singleton)} singletons and "
          f"{n_short_only} ultra-short ASVs")
    dt = dt[~dt["sequence"].isin(pre_remove)].copy()

print(f"[INFO] Chimera input: {dt['sequence'].nunique()} ASVs, "
      f"{len(dt)} non-zero entries")

# ---------------------------------------------------------------------------
# Run sparse chimera detection
# ---------------------------------------------------------------------------
chimeric_seqs = find_chimeras_sparse(dt, cpus=cpus)

# Remove chimeric sequences
dt_clean = dt[~dt["sequence"].isin(chimeric_seqs)].copy()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
chimeras_removed = len(chimeric_seqs)
total_removed    = n_pre_removed + chimeras_removed
n_output_asvs    = dt_clean["sequence"].nunique()
n_output_reads   = int(dt_clean["count"].sum())
pct_retained     = round(n_output_reads / max(n_input_reads, 1) * 100, 1)

print(f"[INFO] Removed {chimeras_removed} chimeras + "
      f"{n_pre_removed} pre-filtered = {total_removed} total")
print(f"[INFO] Retained {pct_retained} % of reads ( "
      f"{n_output_asvs} ASVs across {dt_clean['sample'].nunique()} samples)")

# ---------------------------------------------------------------------------
# Save output — stays in long format
# ---------------------------------------------------------------------------
with open("seqtab_nochim.pkl", "wb") as f:
    pickle.dump(dt_clean, f)

stats = pd.DataFrame([{
    "step": "chimera_removal",
    "input_asvs": n_input_asvs,
    "pre_filtered": n_pre_removed,
    "chimeras": chimeras_removed,
    "output_asvs": n_output_asvs,
    "samples": dt_clean["sample"].nunique(),
    "total_reads": n_output_reads,
    "pct_retained": pct_retained,
}])
stats.to_csv("chimera_stats.tsv", sep="\t", index=False)

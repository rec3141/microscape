#!/usr/bin/env python3
#
# dada2_denoise.py — Denoise, merge pairs, and build a sequence table (per plate)
#
# Python mirror of dada2_denoise.R. For each sample on a plate:
#   1. Dereplicate forward and reverse filtered reads via dada2gpu
#   2. Apply the DADA2 denoising algorithm using pre-learned error models
#   3. Merge denoised forward/reverse pairs (overlap alignment)
#   4. Combine all merged samples into a sequence table (pandas DataFrame)
#   5. Per-plate chimera removal (abundance-ratio bimera check)
#
# Usage:
#   dada2_denoise.py <plate_id> <errF.pkl> <errR.pkl> <min_overlap> <cpus>
#
# Inputs (discovered automatically):
#   *_R1.filt.fastq.gz   Forward filtered reads
#   *_R2.filt.fastq.gz   Reverse filtered reads
#
# Outputs:
#   <plate_id>.seqtab.pkl   Sequence table (pickle, long-format DataFrame)
#   <plate_id>.seqtab.tsv   Sequence table (tab-delimited, for inspection)

import sys
import os
import re
import glob
import pickle
import numpy as np
import pandas as pd

sys.path.insert(0, "/data/dada2_gpu")
import dada2gpu


# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if len(sys.argv) < 6:
    print("Usage: dada2_denoise.py <plate_id> <errF.pkl> <errR.pkl> "
          "<min_overlap> <cpus>", file=sys.stderr)
    sys.exit(1)

plate_id    = sys.argv[1]
errF_path   = sys.argv[2]
errR_path   = sys.argv[3]
min_overlap = int(sys.argv[4])
cpus        = int(sys.argv[5])

# ---------------------------------------------------------------------------
# Load pre-learned error models
# ---------------------------------------------------------------------------
with open(errF_path, "rb") as f:
    errF = pickle.load(f)
with open(errR_path, "rb") as f:
    errR = pickle.load(f)

# ---------------------------------------------------------------------------
# Discover filtered FASTQ files
# ---------------------------------------------------------------------------
fwd_files = sorted(glob.glob("*_R1.filt.fastq.gz"))
rev_files = sorted(glob.glob("*_R2.filt.fastq.gz"))

# Remove empty/missing files (keep pairs aligned)
valid_pairs = []
for f, r in zip(fwd_files, rev_files):
    f_ok = os.path.exists(f) and os.path.getsize(f) > 100
    r_ok = os.path.exists(r) and os.path.getsize(r) > 100
    if f_ok and r_ok:
        valid_pairs.append((f, r))

fwd_files = [p[0] for p in valid_pairs]
rev_files = [p[1] for p in valid_pairs]

# Derive sample names by stripping the _R1.filt.fastq.gz suffix
sample_names = [re.sub(r"_R1\.filt\.fastq\.gz$", "", os.path.basename(f))
                for f in fwd_files]

if len(fwd_files) == 0:
    print("[ERROR] No valid filtered files found for denoising", file=sys.stderr)
    sys.exit(1)

print(f"[INFO] Plate {plate_id} : denoising {len(fwd_files)} samples")


# ---------------------------------------------------------------------------
# Merge pairs implementation
# ---------------------------------------------------------------------------
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def find_best_overlap(fwd_seq, rev_rc, min_ov):
    """Find the best overlap between the 3' end of fwd and 5' end of rev_rc.

    Slides rev_rc along fwd to find the position with the most matches.
    Returns (overlap_length, n_mismatches, merged_seq) or None if no
    valid overlap is found.
    """
    fwd_len = len(fwd_seq)
    rev_len = len(rev_rc)
    best_overlap = 0
    best_mismatches = float('inf')
    best_pos = -1

    # The overlap region: the end of fwd overlaps with the start of rev_rc.
    # overlap_len ranges from min_ov to min(fwd_len, rev_len).
    max_overlap = min(fwd_len, rev_len)

    for ov in range(min_ov, max_overlap + 1):
        fwd_region = fwd_seq[fwd_len - ov:]
        rev_region = rev_rc[:ov]

        mismatches = sum(1 for a, b in zip(fwd_region, rev_region) if a != b)
        match_ratio = 1.0 - mismatches / ov

        if match_ratio > 0.9 and (ov > best_overlap or
                                   (ov == best_overlap and
                                    mismatches < best_mismatches)):
            best_overlap = ov
            best_mismatches = mismatches
            best_pos = ov

    if best_pos < 0:
        return None

    # Build merged sequence: fwd non-overlap + overlap (use fwd) + rev non-overlap
    ov = best_overlap
    merged = fwd_seq + rev_rc[ov:]
    return (ov, best_mismatches, merged)


def merge_pairs(fwd_denoised, rev_denoised, min_ov=10):
    """Merge denoised forward and reverse reads by overlap alignment.

    Args:
        fwd_denoised: dict {sequence: abundance} from dada() on forward reads
        rev_denoised: dict {sequence: abundance} from dada() on reverse reads
        min_ov: minimum overlap length

    Returns:
        dict {merged_sequence: abundance}
    """
    merged = {}

    fwd_seqs = sorted(fwd_denoised.keys(),
                       key=lambda s: fwd_denoised[s], reverse=True)
    rev_seqs = sorted(rev_denoised.keys(),
                       key=lambda s: rev_denoised[s], reverse=True)

    for fwd_seq in fwd_seqs:
        fwd_abund = fwd_denoised[fwd_seq]
        for rev_seq in rev_seqs:
            rev_abund = rev_denoised[rev_seq]

            # Use the minimum abundance of the pair
            pair_abund = min(fwd_abund, rev_abund)
            if pair_abund <= 0:
                continue

            rev_rc = reverse_complement(rev_seq)
            result = find_best_overlap(fwd_seq, rev_rc, min_ov)

            if result is not None:
                ov, mm, merged_seq = result
                merged[merged_seq] = merged.get(merged_seq, 0) + pair_abund

    return merged


# ---------------------------------------------------------------------------
# Per-plate chimera removal (bimera detection)
#
# A simple abundance-ratio approach: for each ASV, check if it can be
# reconstructed as a chimera of two more-abundant parent ASVs. A sequence
# is flagged as a bimera if its left half closely matches one parent and
# its right half closely matches another.
# ---------------------------------------------------------------------------
def is_bimera(query, parents, min_fold=1.5, min_match_frac=0.9):
    """Check if query is a bimera (chimera of two parents).

    Splits the query at every possible breakpoint and checks if the left
    half matches one parent and the right half matches another parent.
    """
    qlen = len(query)
    if qlen < 20:
        return False

    for bp in range(10, qlen - 10):
        left = query[:bp]
        right = query[bp:]
        left_len = len(left)
        right_len = len(right)

        for i, p1 in enumerate(parents):
            if len(p1) < qlen:
                continue
            # Check left half against p1
            p1_left = p1[:left_len]
            left_matches = sum(1 for a, b in zip(left, p1_left) if a == b)
            if left_matches / left_len < min_match_frac:
                continue

            for j, p2 in enumerate(parents):
                if i == j:
                    continue
                if len(p2) < qlen:
                    continue
                # Check right half against p2 (aligned from the right end)
                p2_right = p2[len(p2) - right_len:]
                right_matches = sum(1 for a, b in zip(right, p2_right)
                                    if a == b)
                if right_matches / right_len >= min_match_frac:
                    return True

    return False


def remove_bimeras_from_table(seqtab_dict, min_fold=1.5):
    """Remove bimeric sequences from a sequence table.

    Args:
        seqtab_dict: dict {sequence: total_abundance}

    Returns:
        set of chimeric sequences
    """
    # Sort by abundance descending
    sorted_seqs = sorted(seqtab_dict.keys(),
                          key=lambda s: seqtab_dict[s], reverse=True)

    chimeras = set()

    for idx, query in enumerate(sorted_seqs):
        query_abund = seqtab_dict[query]
        # Parents must be more abundant (by min_fold)
        parents = [s for s in sorted_seqs[:idx]
                   if seqtab_dict[s] >= min_fold * query_abund]

        if len(parents) < 2:
            continue

        # Limit to top 20 parents to keep runtime reasonable
        parents = parents[:20]

        if is_bimera(query, parents, min_fold=min_fold):
            chimeras.add(query)

    return chimeras


# ---------------------------------------------------------------------------
# Per-sample: dereplicate -> denoise -> merge pairs
# ---------------------------------------------------------------------------
all_merged = {}  # sample_name -> {merged_seq: abundance}

for i, sample_name in enumerate(sample_names):
    print(f"[INFO] Processing: {sample_name}")

    # Dereplicate: collapse identical reads, track abundances
    derepF = dada2gpu.derep_fastq(fwd_files[i])
    derepR = dada2gpu.derep_fastq(rev_files[i])

    # Denoise: apply error model to distinguish real variants from errors
    try:
        ddF = dada2gpu.dada(derepF, err=errF)
        ddR = dada2gpu.dada(derepR, err=errR)
    except Exception as e:
        print(f"[WARNING] {sample_name}: dada() failed: {e}")
        continue

    # Extract denoised sequences (may be missing if sample is too small)
    fwd_denoised = ddF.get("denoised", {})
    rev_denoised = ddR.get("denoised", {})

    if not fwd_denoised or not rev_denoised:
        print(f"[WARNING] {sample_name}: no denoised sequences, skipping")
        continue

    merged = merge_pairs(fwd_denoised, rev_denoised, min_ov=min_overlap)

    if len(merged) > 0:
        all_merged[sample_name] = merged

# Drop samples that produced no merged reads
if len(all_merged) == 0:
    print(f"[ERROR] No merged reads produced for plate {plate_id}",
          file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Build sequence table
# ---------------------------------------------------------------------------
# Collect all data in long format
rows = []
for sample_name, seq_dict in all_merged.items():
    for seq, count in seq_dict.items():
        rows.append({"sample": sample_name, "sequence": seq, "count": int(count)})

dt = pd.DataFrame(rows)
# Aggregate any duplicates
dt = dt.groupby(["sample", "sequence"], as_index=False)["count"].sum()

n_samples = dt["sample"].nunique()
n_asvs = dt["sequence"].nunique()
total_reads = dt["count"].sum()

print(f"[INFO] Plate {plate_id} : raw table has {n_samples} samples, "
      f"{n_asvs} unique ASVs, {total_reads} total reads")

# ---------------------------------------------------------------------------
# Per-plate chimera removal
# ---------------------------------------------------------------------------
seq_totals = dt.groupby("sequence")["count"].sum().to_dict()

chimeric_seqs = remove_bimeras_from_table(seq_totals)

dt_nochim = dt[~dt["sequence"].isin(chimeric_seqs)].copy()

n_chimeras = len(chimeric_seqs)
n_asvs_after = dt_nochim["sequence"].nunique()
reads_after = dt_nochim["count"].sum()
pct = round(reads_after / max(total_reads, 1) * 100, 1)

print(f"[INFO] Plate {plate_id} : removed {n_chimeras} chimeras, "
      f"retained {pct} % of reads ( {n_asvs_after} ASVs)")

# ---------------------------------------------------------------------------
# Save in both formats: pickle for downstream Python steps, TSV for inspection
# ---------------------------------------------------------------------------
with open(f"{plate_id}.seqtab.pkl", "wb") as f:
    pickle.dump(dt_nochim, f)

# Wide TSV for inspection
wide = dt_nochim.pivot_table(index="sample", columns="sequence",
                              values="count", fill_value=0, aggfunc="sum")
wide.to_csv(f"{plate_id}.seqtab.tsv", sep="\t")

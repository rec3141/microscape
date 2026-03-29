#!/usr/bin/env python3
#
# dada2_denoise.py — Denoise, merge pairs, and build a sequence table (per plate)
#
# Python mirror of dada2_denoise.R. For each sample on a plate:
#   1. Dereplicate forward and reverse filtered reads via dada2py
#   2. Apply the DADA2 denoising algorithm using pre-learned error models
#   3. Merge denoised forward/reverse pairs via C-level NW alignment
#   4. Combine all merged samples into a sequence table (pandas DataFrame)
#   5. Per-plate chimera removal via C-level consensus bimera detection
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

import py as dada2py


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
# Merge pairs and chimera removal use dada2's C-level implementations
# via the dada2py library (py/paired.py and py/chimera.py).
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Per-sample: dereplicate -> denoise -> merge pairs
# ---------------------------------------------------------------------------
all_merged = {}  # sample_name -> {merged_seq: abundance}

for i, sample_name in enumerate(sample_names):
    print(f"[INFO] Processing: {sample_name}")

    # Dereplicate: collapse identical reads, track abundances
    derepF = dada2py.derep_fastq(fwd_files[i])
    derepR = dada2py.derep_fastq(rev_files[i])

    # Denoise: apply error model to distinguish real variants from errors
    try:
        ddF = dada2py.dada(derepF, err=errF)
        ddR = dada2py.dada(derepR, err=errR)
    except Exception as e:
        print(f"[WARNING] {sample_name}: dada() failed: {e}")
        continue

    # Check for denoised output
    if not ddF.get("cluster_seqs") or not ddR.get("cluster_seqs"):
        print(f"[WARNING] {sample_name}: no denoised sequences, skipping")
        continue

    # Merge pairs using C-level NW alignment (matches R's mergePairs)
    merged_list = dada2py.merge_pairs(ddF, derepF, ddR, derepR,
                                       min_overlap=min_overlap,
                                       trim_overhang=True, verbose=False)

    # Collect accepted merged reads into {sequence: abundance}
    accepted = {}
    for m in merged_list:
        if m["accept"]:
            accepted[m["sequence"]] = accepted.get(m["sequence"], 0) + m["abundance"]

    if accepted:
        all_merged[sample_name] = accepted

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
# Per-plate chimera removal using dada2's C-level consensus method
# ---------------------------------------------------------------------------
all_seqs = sorted(set(s for d in all_merged.values() for s in d))
sample_list = list(all_merged.keys())
seq_to_idx = {s: i for i, s in enumerate(all_seqs)}

mat = np.zeros((len(sample_list), len(all_seqs)), dtype=np.int32)
for si, sam in enumerate(sample_list):
    for seq, ab in all_merged[sam].items():
        mat[si, seq_to_idx[seq]] = ab

seqtab = {"table": mat, "seqs": all_seqs}
chim_result = dada2py.remove_bimera_denovo(seqtab, verbose=True)

is_chimera = np.array(chim_result["is_chimera"])
chimeric_seqs = set(s for s, c in zip(all_seqs, is_chimera) if c)
dt_nochim = dt[~dt["sequence"].isin(chimeric_seqs)].copy()

n_chimeras = int(is_chimera.sum())
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

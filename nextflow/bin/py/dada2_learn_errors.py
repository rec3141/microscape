#!/usr/bin/env python3
#
# dada2_learn_errors.py — Learn error rates from filtered reads (per plate)
#
# Python mirror of dada2_learn_errors.R. Uses dada2gpu.learn_errors() to
# learn per-plate error models from filtered FASTQ files. Saves error
# matrices as pickle files and generates diagnostic error-rate plots.
#
# Usage:
#   dada2_learn_errors.py <plate_id> <cpus>
#
# Inputs (discovered automatically):
#   *_R1.filt.fastq.gz   Forward filtered reads (one per sample)
#   *_R2.filt.fastq.gz   Reverse filtered reads (one per sample)
#
# Outputs:
#   <plate_id>_errF.pkl          Forward error model (numpy array)
#   <plate_id>_errR.pkl          Reverse error model (numpy array)
#   <plate_id>_error_rates.pdf   Diagnostic plots (observed vs fitted rates)

import sys
import os
import glob
import pickle
import numpy as np

sys.path.insert(0, "/data/dada2_gpu")
import dada2gpu

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if len(sys.argv) < 3:
    print("Usage: dada2_learn_errors.py <plate_id> <cpus>", file=sys.stderr)
    sys.exit(1)

plate_id = sys.argv[1]
cpus     = int(sys.argv[2])

# ---------------------------------------------------------------------------
# Discover filtered FASTQ files in the working directory
# ---------------------------------------------------------------------------
fwd_files = sorted(glob.glob("*_R1.filt.fastq.gz"))
rev_files = sorted(glob.glob("*_R2.filt.fastq.gz"))

# Drop any empty or missing files (can happen when all reads are filtered out)
def filter_valid(files, min_size=100):
    """Keep only files that exist and are larger than min_size bytes.
    A gzip file with zero reads is still ~20 bytes (header only), so
    we use 100 as threshold to skip truly empty files."""
    valid = []
    for f in files:
        try:
            sz = os.path.getsize(f)
            if sz > min_size:
                valid.append(f)
            else:
                print(f"[WARNING] Skipping empty file: {f} ({sz} bytes)")
        except OSError:
            pass
    return valid

fwd_files = filter_valid(fwd_files)
rev_files = filter_valid(rev_files)

if len(fwd_files) == 0:
    print(f"[ERROR] No valid filtered FASTQ files found for plate {plate_id}",
          file=sys.stderr)
    sys.exit(1)

print(f"[INFO] Plate {plate_id} : learning error rates from "
      f"{len(fwd_files)} samples")

# ---------------------------------------------------------------------------
# Learn error models
#
# dada2gpu.learn_errors() pools quality-score information across all input
# files and fits a parametric error model. Returns a 16 x ncol numpy array
# representing transition probabilities for each quality score.
# ---------------------------------------------------------------------------
errF = dada2gpu.learn_errors(fwd_files, nbases=1e8)
errR = dada2gpu.learn_errors(rev_files, nbases=1e8)

with open(f"{plate_id}_errF.pkl", "wb") as f:
    pickle.dump(errF, f)
with open(f"{plate_id}_errR.pkl", "wb") as f:
    pickle.dump(errR, f)

# ---------------------------------------------------------------------------
# Diagnostic plots — observed error rates vs fitted model
#
# The error matrix is 16 x ncol where rows correspond to the 16 possible
# nucleotide transitions (A2A, A2C, A2G, A2T, C2A, ..., T2T) and columns
# correspond to quality scores.
# ---------------------------------------------------------------------------
TRANSITIONS = [
    "A→A", "A→C", "A→G", "A→T",
    "C→A", "C→C", "C→G", "C→T",
    "G→A", "G→C", "G→G", "G→T",
    "T→A", "T→C", "T→G", "T→T",
]

def plot_error_model(err_matrix, title, ax_array):
    """Plot the 16-panel error model on a 4x4 grid of axes."""
    nrows, ncols = err_matrix.shape
    qual_scores = np.arange(ncols)

    for idx in range(min(16, nrows)):
        row = idx // 4
        col = idx % 4
        ax = ax_array[row, col]
        rates = err_matrix[idx, :]
        # Plot log10 error rate vs quality score
        with np.errstate(divide='ignore'):
            log_rates = np.log10(np.maximum(rates, 1e-10))
        ax.plot(qual_scores, log_rates, 'o-', markersize=2, linewidth=0.8,
                color='steelblue')
        # Nominal error rate line
        nominal = np.log10(np.maximum(10.0 ** (-qual_scores / 10.0), 1e-10))
        ax.plot(qual_scores, nominal, '--', color='red', linewidth=0.5,
                alpha=0.5)
        ax.set_title(TRANSITIONS[idx], fontsize=8)
        ax.set_ylim(-8, 0)
        ax.tick_params(labelsize=6)

    ax_array[3, 0].set_xlabel("Quality Score", fontsize=8)
    ax_array[3, 0].set_ylabel("log10(Error Rate)", fontsize=8)


fig, axes = plt.subplots(4, 4, figsize=(10, 10))
fig.suptitle(f"{plate_id} — Forward Error Model", fontsize=12)
plot_error_model(errF, "Forward", axes)
fig.tight_layout(rect=[0, 0, 1, 0.96])

fig2, axes2 = plt.subplots(4, 4, figsize=(10, 10))
fig2.suptitle(f"{plate_id} — Reverse Error Model", fontsize=12)
plot_error_model(errR, "Reverse", axes2)
fig2.tight_layout(rect=[0, 0, 1, 0.96])

from matplotlib.backends.backend_pdf import PdfPages
with PdfPages(f"{plate_id}_error_rates.pdf") as pdf:
    pdf.savefig(fig)
    pdf.savefig(fig2)

plt.close('all')

print(f"[INFO] Plate {plate_id} : error models saved")

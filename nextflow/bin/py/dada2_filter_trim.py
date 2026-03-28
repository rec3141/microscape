#!/usr/bin/env python3
#
# dada2_filter_trim.py — Per-sample quality filtering of paired-end reads
#
# Python mirror of dada2_filter_trim.R. Uses BioPython for FASTQ parsing.
# Removes low-quality reads, truncates at first base <= truncQ, optionally
# truncates to fixed lengths, and filters by expected errors and N count.
#
# Usage:
#   dada2_filter_trim.py <sample_id> <fwd.fastq.gz> <rev.fastq.gz> \
#       <maxEE> <truncQ> <maxN> <truncLen_fwd> <truncLen_rev> <cpus>
#
# Outputs:
#   <sample_id>_R1.filt.fastq.gz   Filtered forward reads
#   <sample_id>_R2.filt.fastq.gz   Filtered reverse reads
#   <sample_id>_filt_stats.tsv      Filtering statistics

import sys
import os
import gzip
import math
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
if len(sys.argv) < 10:
    print("Usage: dada2_filter_trim.py <sample_id> <r1> <r2> "
          "<maxEE> <truncQ> <maxN> <truncLen_fwd> <truncLen_rev> <cpus>",
          file=sys.stderr)
    sys.exit(1)

sample_id    = sys.argv[1]
r1_in        = sys.argv[2]
r2_in        = sys.argv[3]
max_ee       = float(sys.argv[4])
trunc_q      = int(sys.argv[5])
max_n        = int(sys.argv[6])
trunc_fwd    = int(sys.argv[7])
trunc_rev    = int(sys.argv[8])
cpus         = int(sys.argv[9])

# Output filenames
r1_out = f"{sample_id}_R1.filt.fastq.gz"
r2_out = f"{sample_id}_R2.filt.fastq.gz"


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def expected_errors(quality_scores):
    """Calculate expected errors = sum(10^(-Q/10)) for a read."""
    return sum(10.0 ** (-q / 10.0) for q in quality_scores)


def count_ns(seq_str):
    """Count the number of N bases in a sequence."""
    return seq_str.upper().count('N')


def truncate_at_quality(record, min_qual):
    """Truncate read at the first base with quality <= min_qual.

    Returns the truncated SeqRecord, or the original if no base falls
    at or below the threshold.
    """
    quals = record.letter_annotations["phred_quality"]
    for i, q in enumerate(quals):
        if q <= min_qual:
            return record[:i]
    return record


def filter_read(record, max_n_val, max_ee_val, trunc_q_val, trunc_len):
    """Apply the filter cascade to a single read.

    Returns the filtered SeqRecord, or None if the read fails.
    Steps:
      1. Truncate at first base with quality <= trunc_q
      2. Optionally truncate to fixed length
      3. Reject if too many Ns
      4. Reject if expected errors > maxEE
      5. Reject if empty after truncation
    """
    # Truncate at low quality
    rec = truncate_at_quality(record, trunc_q_val)

    # Fixed-length truncation
    if trunc_len > 0:
        if len(rec) < trunc_len:
            return None  # too short after quality truncation
        rec = rec[:trunc_len]

    # Reject reads too short for kmer-based processing
    # R's filterAndTrim uses minLen=20 by default
    if len(rec) < 20:
        return None

    # Reject reads with too many Ns
    if count_ns(str(rec.seq)) > max_n_val:
        return None

    # Reject reads with too many expected errors
    quals = rec.letter_annotations["phred_quality"]
    if expected_errors(quals) > max_ee_val:
        return None

    return rec


# ---------------------------------------------------------------------------
# Filter and trim
# ---------------------------------------------------------------------------
trunc_len_str = f"{trunc_fwd},{trunc_rev}"
print(f"[INFO] {sample_id} : filtering with maxEE={max_ee} "
      f"truncQ={trunc_q} truncLen={trunc_len_str}")

reads_in = 0
reads_out = 0

# Open input files (gzipped FASTQ)
fwd_handle = gzip.open(r1_in, "rt")
rev_handle = gzip.open(r2_in, "rt")

# Open output files (gzipped FASTQ)
fwd_out_handle = gzip.open(r1_out, "wt")
rev_out_handle = gzip.open(r2_out, "wt")

try:
    fwd_iter = SeqIO.parse(fwd_handle, "fastq")
    rev_iter = SeqIO.parse(rev_handle, "fastq")

    for fwd_rec, rev_rec in zip(fwd_iter, rev_iter):
        reads_in += 1

        fwd_filt = filter_read(fwd_rec, max_n, max_ee, trunc_q, trunc_fwd)
        rev_filt = filter_read(rev_rec, max_n, max_ee, trunc_q, trunc_rev)

        # Both reads must pass for the pair to be kept
        if fwd_filt is not None and rev_filt is not None:
            SeqIO.write(fwd_filt, fwd_out_handle, "fastq")
            SeqIO.write(rev_filt, rev_out_handle, "fastq")
            reads_out += 1

finally:
    fwd_handle.close()
    rev_handle.close()
    fwd_out_handle.close()
    rev_out_handle.close()

# Handle samples where all reads were filtered out
# Match R's filterAndTrim behavior: remove output files and recreate as 0-byte
# so downstream steps can detect by file size (>20 bytes = has reads)
if reads_out == 0:
    print(f"[WARNING] {sample_id} : no reads passed filter")
    os.remove(r1_out)
    os.remove(r2_out)
    # Create 0-byte files (Nextflow requires declared outputs to exist)
    open(r1_out, 'w').close()
    open(r2_out, 'w').close()

pct = round(reads_out / max(reads_in, 1) * 100, 1)

# ---------------------------------------------------------------------------
# Write stats TSV
# ---------------------------------------------------------------------------
stats_path = f"{sample_id}_filt_stats.tsv"
with open(stats_path, "w") as f:
    f.write("sample\treads_in\treads_out\tpct_retained\n")
    f.write(f"{sample_id}\t{reads_in}\t{reads_out}\t{pct}\n")

print(f"[INFO] {sample_id} : {reads_out} / {reads_in} "
      f"reads passed filter ( {pct} %)")

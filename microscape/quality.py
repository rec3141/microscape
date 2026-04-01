"""Quality profiling and automated truncation length selection.

Analyzes per-position quality scores from paired-end FASTQ files and
recommends truncation lengths for downstream filtering.
"""

import gzip
import glob
import os
from typing import List, Optional, Tuple

import numpy as np


def _read_quals(files: List[str], n_reads: int) -> List[List[int]]:
    """Read quality scores from FASTQ files."""
    all_quals = []
    count = 0
    for fpath in files:
        opener = gzip.open if fpath.endswith(".gz") else open
        with opener(fpath, "rt") as fh:
            while count < n_reads:
                header = fh.readline()
                if not header:
                    break
                fh.readline()  # seq
                fh.readline()  # +
                qual_str = fh.readline().rstrip()
                if not qual_str:
                    break
                all_quals.append([ord(c) - 33 for c in qual_str])
                count += 1
        if count >= n_reads:
            break
    return all_quals


def _find_trunc_pos(quals_list: List[List[int]], min_q: float,
                    window: int) -> int:
    """Find position where rolling median quality drops below threshold."""
    if not quals_list:
        return 0
    max_len = max(len(q) for q in quals_list)
    mat = np.full((len(quals_list), max_len), np.nan)
    for i, q in enumerate(quals_list):
        mat[i, :len(q)] = q
    medians = np.nanmedian(mat, axis=0)
    if len(medians) < window:
        return int(len(medians))
    rolling = np.convolve(medians, np.ones(window) / window, mode="valid")
    for i, val in enumerate(rolling):
        if val < min_q:
            return i + window // 2
    return int(len(medians))


def auto_trim(
    input_dir: str,
    *,
    min_quality: float = 25.0,
    window: int = 10,
    n_reads: int = 10000,
    n_files: int = 20,
    verbose: bool = False,
) -> dict:
    """Analyze quality profiles and recommend truncation lengths.

    Samples reads from paired-end FASTQs in *input_dir*, computes
    per-position median quality, and finds where a rolling median drops
    below *min_quality*.

    Args:
        input_dir: Directory containing paired FASTQ files.
        min_quality: Minimum rolling median Q score to keep a position.
        window: Rolling window size for quality smoothing.
        n_reads: Number of reads to sample.
        n_files: Maximum number of file pairs to sample from.
        verbose: Print progress messages.

    Returns:
        Dict with keys: trunc_len_fwd, trunc_len_rev, fwd_read_len,
        rev_read_len, n_reads_sampled, min_quality.
    """
    fwd_files = sorted(
        glob.glob(os.path.join(input_dir, "*_1.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*_R1*.fastq.gz"))
    )
    rev_files = sorted(
        glob.glob(os.path.join(input_dir, "*_2.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*_R2*.fastq.gz"))
    )

    if not fwd_files or not rev_files:
        raise FileNotFoundError(
            f"No paired FASTQ files found in {input_dir}"
        )

    # Sample from a subset of files
    n_sample = min(n_files, len(fwd_files))
    step = max(1, len(fwd_files) // n_sample)
    fwd_sample = fwd_files[::step][:n_sample]
    rev_sample = rev_files[::step][:n_sample]

    if verbose:
        print(f"[INFO] Sampling {n_reads} reads from {n_sample} file pairs")

    fwd_quals = _read_quals(fwd_sample, n_reads)
    rev_quals = _read_quals(rev_sample, n_reads)

    trunc_fwd = _find_trunc_pos(fwd_quals, min_quality, window)
    trunc_rev = _find_trunc_pos(rev_quals, min_quality, window)

    fwd_len = max(len(q) for q in fwd_quals) if fwd_quals else 0
    rev_len = max(len(q) for q in rev_quals) if rev_quals else 0

    if verbose:
        fwd_end_q = float(np.median([q[-1] for q in fwd_quals if q])) if fwd_quals else 0
        rev_end_q = float(np.median([q[-1] for q in rev_quals if q])) if rev_quals else 0
        print(f"[INFO] Forward: len={fwd_len}, Q@end={fwd_end_q:.0f}, trunc={trunc_fwd}")
        print(f"[INFO] Reverse: len={rev_len}, Q@end={rev_end_q:.0f}, trunc={trunc_rev}")

    return {
        "trunc_len_fwd": trunc_fwd,
        "trunc_len_rev": trunc_rev,
        "fwd_read_len": fwd_len,
        "rev_read_len": rev_len,
        "n_reads_sampled": len(fwd_quals),
        "min_quality": min_quality,
    }

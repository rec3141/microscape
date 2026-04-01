"""Prepare reads for pipeline input.

Reads a metadata file with optional ``dada_run`` and ``dada_plate`` columns,
matches samples to FASTQ files, and creates symlinks with standardized
names that encode run/plate grouping for per-plate error learning.

Output filenames: ``{run}_{plate}_{hash}_R{1,2}.fastq.gz``
where hash is a short unique identifier derived from the original sample ID.
"""

import hashlib
import os
import re
from pathlib import Path
from typing import Optional

import pandas as pd


def _clean_id(s: str) -> str:
    """Sanitize a string for use in filenames."""
    return re.sub(r"[^A-Za-z0-9._-]", "", str(s).strip())


def _short_hash(s: str, length: int = 6) -> str:
    """Deterministic short hash of a string."""
    return hashlib.md5(s.encode()).hexdigest()[:length]


def _find_fastq_pair(reads_dir: Path, sample_id: str):
    """Find R1/R2 or _1/_2 FASTQ files matching a sample ID."""
    patterns = [
        (f"{sample_id}_R1*.fastq.gz", f"{sample_id}_R2*.fastq.gz"),
        (f"{sample_id}_1.fastq.gz", f"{sample_id}_2.fastq.gz"),
        (f"{sample_id}_R1*.fq.gz", f"{sample_id}_R2*.fq.gz"),
    ]
    for p1, p2 in patterns:
        r1 = sorted(reads_dir.glob(p1))
        r2 = sorted(reads_dir.glob(p2))
        if r1 and r2:
            return r1[0], r2[0]
    return None, None


def prep_reads(
    metadata_path: str,
    reads_dir: str,
    output_dir: str,
    *,
    sample_col: str = "Run",
    run_col: str = "dada_run",
    plate_col: str = "dada_plate",
    default_run: str = "run1",
    default_plate: str = "plate1",
    verbose: bool = False,
) -> pd.DataFrame:
    """Create symlinks with standardized names encoding run/plate grouping.

    Args:
        metadata_path: Path to metadata CSV/TSV.
        reads_dir: Directory containing input FASTQ files.
        output_dir: Directory for symlinks.
        sample_col: Column matching sample IDs to FASTQ filenames.
        run_col: Column for sequencing run grouping. If absent, uses
            *default_run* for all samples.
        plate_col: Column for plate grouping within a run. If absent,
            uses *default_plate* for all samples.
        default_run: Default run name when *run_col* is missing.
        default_plate: Default plate name when *plate_col* is missing.
        verbose: Print progress.

    Returns:
        DataFrame with columns: sample_id, run, plate, group, hash,
        new_name, r1_orig, r2_orig.
    """
    reads_path = Path(reads_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Read metadata — detect separator
    if metadata_path.endswith(".tsv") or metadata_path.endswith(".txt"):
        meta = pd.read_csv(metadata_path, sep="\t")
    else:
        meta = pd.read_csv(metadata_path)

    if sample_col not in meta.columns:
        raise ValueError(
            f"Sample column '{sample_col}' not found in metadata. "
            f"Available: {list(meta.columns)}"
        )

    has_run = run_col in meta.columns
    has_plate = plate_col in meta.columns

    if verbose:
        print(f"[INFO] {len(meta)} samples in metadata")
        print(f"[INFO] run column '{run_col}': {'found' if has_run else 'not found, using ' + repr(default_run)}")
        print(f"[INFO] plate column '{plate_col}': {'found' if has_plate else 'not found, using ' + repr(default_plate)}")

    records = []
    linked = 0
    skipped = 0

    for _, row in meta.iterrows():
        sample_id = str(row[sample_col]).strip()
        run = _clean_id(row[run_col]) if has_run else default_run
        plate = _clean_id(row[plate_col]) if has_plate else default_plate
        group = f"{run}_{plate}"
        h = _short_hash(sample_id)

        r1_orig, r2_orig = _find_fastq_pair(reads_path, sample_id)
        if r1_orig is None:
            if verbose:
                print(f"[WARN] No FASTQ files found for {sample_id}")
            skipped += 1
            continue

        new_base = f"{group}_{h}"
        r1_link = out_path / f"{new_base}_R1.fastq.gz"
        r2_link = out_path / f"{new_base}_R2.fastq.gz"

        # Create symlinks (overwrite existing)
        for link, target in [(r1_link, r1_orig), (r2_link, r2_orig)]:
            if link.exists() or link.is_symlink():
                link.unlink()
            link.symlink_to(target.resolve())

        records.append({
            "sample_id": sample_id,
            "run": run,
            "plate": plate,
            "group": group,
            "hash": h,
            "new_name": new_base,
            "r1_orig": str(r1_orig),
            "r2_orig": str(r2_orig),
        })
        linked += 1

    mapping = pd.DataFrame(records)

    # Write mapping TSV
    mapping_path = out_path / "sample_mapping.tsv"
    mapping.to_csv(mapping_path, sep="\t", index=False)

    if verbose:
        groups = mapping["group"].nunique() if len(mapping) > 0 else 0
        print(f"[INFO] Linked {linked} samples, skipped {skipped}")
        print(f"[INFO] {groups} unique run_plate groups")
        print(f"[INFO] Mapping saved to {mapping_path}")

    return mapping

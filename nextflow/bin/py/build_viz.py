#!/usr/bin/env python3
"""
build_viz.py - Convert pipeline outputs (pickle files) into JSON files for the
Svelte frontend visualization.

Replaces the R build_shiny.R script.

Usage:
  build_viz.py <seqtab.pkl> <renorm.pkl> <taxonomy_dir> <metadata.pkl_or_NONE> \
               <sample_tsne.pkl> <seq_tsne.pkl> <network.pkl>

Outputs (written to current working directory):
  samples.json        Sample metadata with t-SNE coordinates
  asvs.json.gz        ASV info with t-SNE coordinates and taxonomy
  counts.json.gz      Sparse count matrix (gzipped)
  network.json        Correlation edge list
  taxonomy.json       Per-database taxonomy assignments
  renorm_stats.json   Group-level summary statistics
"""

import sys
import os
import json
import gzip
import glob
import pickle
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def log_info(msg):
    print(f"[INFO] {msg}", flush=True)


def log_warn(msg):
    print(f"[WARN] {msg}", flush=True)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_pickle(path):
    """Load a pickle file and return its contents."""
    log_info(f"Loading {path}")
    with open(path, "rb") as fh:
        return pickle.load(fh)


def write_json(obj, path, compress=False):
    """Write an object as JSON, optionally gzip-compressed.

    When compress=True, writes both .json.gz and .json versions.
    When compress=False, writes only .json.
    """
    json_str = json.dumps(obj, separators=(",", ":"), allow_nan=False)

    if compress:
        gz_path = path if path.endswith(".gz") else path + ".gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
            fh.write(json_str)
        log_info(f"Wrote {gz_path} ({os.path.getsize(gz_path):,} bytes gzipped)")
    else:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(json_str)
        log_info(f"Wrote {path} ({os.path.getsize(path):,} bytes)")


def round_float(val, decimals=4):
    """Round a float value, handling NaN/Inf gracefully."""
    if val is None or (isinstance(val, float) and (np.isnan(val) or np.isinf(val))):
        return 0.0
    return round(float(val), decimals)


# ---------------------------------------------------------------------------
# Build samples.json
# ---------------------------------------------------------------------------

def build_samples(seqtab, sample_tsne, metadata):
    """Build the samples array with t-SNE coordinates and metadata.

    Returns:
        list of dicts, one per sample
    """
    log_info("Building samples.json")

    # Aggregate per-sample stats from count table
    sample_stats = seqtab.groupby("sample").agg(
        total_reads=("count", "sum"),
        n_asvs=("sequence", "nunique"),
    ).reset_index()

    # Build t-SNE lookup: label -> (x, y)
    tsne_lookup = {}
    for _, row in sample_tsne.iterrows():
        tsne_lookup[row["label"]] = (
            round_float(row["tSNE1"]),
            round_float(row["tSNE2"]),
        )

    # Build metadata lookup if available
    meta_lookup = {}
    meta_fields = []
    if metadata is not None:
        if isinstance(metadata, pd.DataFrame):
            # Determine the sample ID column
            id_col = None
            for candidate in ("sample", "Sample", "sample_id", "SampleID", "id"):
                if candidate in metadata.columns:
                    id_col = candidate
                    break
            if id_col is None and metadata.index.name:
                # Use index as sample ID
                metadata = metadata.reset_index()
                id_col = metadata.columns[0]
            elif id_col is None:
                id_col = metadata.columns[0]

            meta_fields = [c for c in metadata.columns if c != id_col]
            for _, row in metadata.iterrows():
                sid = str(row[id_col])
                meta_lookup[sid] = {
                    col: _safe_json_value(row[col]) for col in meta_fields
                }
            log_info(f"Metadata: {len(meta_lookup)} samples, fields: {meta_fields}")

    # Assemble sample records
    records = []
    for _, row in sample_stats.iterrows():
        sid = row["sample"]
        x, y = tsne_lookup.get(sid, (0.0, 0.0))
        rec = {
            "id": sid,
            "x": x,
            "y": y,
            "total_reads": int(row["total_reads"]),
            "n_asvs": int(row["n_asvs"]),
        }
        # Merge metadata fields
        if sid in meta_lookup:
            rec.update(meta_lookup[sid])
        records.append(rec)

    log_info(f"samples.json: {len(records)} samples")
    return records


def _safe_json_value(val):
    """Convert a value to a JSON-safe type."""
    if val is None:
        return None
    if isinstance(val, (np.integer,)):
        return int(val)
    if isinstance(val, (np.floating,)):
        if np.isnan(val) or np.isinf(val):
            return None
        return round(float(val), 4)
    if isinstance(val, float):
        if np.isnan(val) or np.isinf(val):
            return None
        return round(val, 4)
    if isinstance(val, (np.bool_,)):
        return bool(val)
    return str(val)


# ---------------------------------------------------------------------------
# Build asvs.json
# ---------------------------------------------------------------------------

def build_asvs(seqtab, seq_tsne, renorm_table_list, taxonomy_dict):
    """Build the ASVs array with t-SNE coordinates, group, and taxonomy.

    Returns:
        tuple of (list of dicts, list of ASV IDs in order)
    """
    log_info("Building asvs.json")

    # Get unique sequences sorted for stable ASV IDs
    sequences = sorted(seqtab["sequence"].unique())
    seq_to_id = {seq: f"ASV_{i+1:04d}" for i, seq in enumerate(sequences)}

    # Aggregate per-ASV stats
    asv_stats = seqtab.groupby("sequence").agg(
        total_reads=("count", "sum"),
        n_samples=("sample", "nunique"),
    ).reset_index()

    # t-SNE lookup
    tsne_lookup = {}
    for _, row in seq_tsne.iterrows():
        tsne_lookup[row["label"]] = (
            round_float(row["tSNE1"]),
            round_float(row["tSNE2"]),
        )

    # Group lookup from renorm
    group_lookup = {}
    if renorm_table_list is not None:
        for _, row in renorm_table_list.iterrows():
            group_lookup[row["sequence"]] = row["group"]

    # Taxonomy lookup: pick first available database for the concatenated string
    tax_string_lookup = {}
    if taxonomy_dict:
        # Use the first database
        first_db = next(iter(taxonomy_dict))
        tax_data = taxonomy_dict[first_db]
        if "tax" in tax_data and isinstance(tax_data["tax"], pd.DataFrame):
            tax_df = tax_data["tax"]
            for seq in tax_df.index:
                vals = [str(v) if pd.notna(v) else "" for v in tax_df.loc[seq]]
                tax_string_lookup[seq] = ";".join(vals)

    # Assemble ASV records
    records = []
    asv_id_list = []
    for _, row in asv_stats.iterrows():
        seq = row["sequence"]
        asv_id = seq_to_id[seq]
        x, y = tsne_lookup.get(seq, (0.0, 0.0))

        rec = {
            "id": asv_id,
            "sequence": seq,
            "x": x,
            "y": y,
            "total_reads": int(row["total_reads"]),
            "n_samples": int(row["n_samples"]),
            "group": group_lookup.get(seq, "unknown"),
            "taxonomy": tax_string_lookup.get(seq, ""),
        }
        records.append(rec)
        asv_id_list.append(asv_id)

    log_info(f"asvs.json: {len(records)} ASVs")
    return records, seq_to_id


# ---------------------------------------------------------------------------
# Build counts.json.gz
# ---------------------------------------------------------------------------

def build_counts(seqtab, seq_to_id):
    """Build the sparse count matrix.

    Returns:
        dict with 'data', 'samples', 'asvs' keys
    """
    log_info("Building counts.json.gz")

    # Ordered sample and ASV lists
    samples = sorted(seqtab["sample"].unique())
    asvs = sorted(seq_to_id.keys(), key=lambda s: seq_to_id[s])
    asv_ids = [seq_to_id[s] for s in asvs]

    sample_idx = {s: i for i, s in enumerate(samples)}
    asv_idx = {s: i for i, s in enumerate(asvs)}

    # Compute per-sample totals for proportions
    sample_totals = seqtab.groupby("sample")["count"].sum()

    # Build sparse data array: [sample_idx, asv_idx, count, proportion]
    data = []
    for _, row in seqtab.iterrows():
        cnt = int(row["count"])
        if cnt == 0:
            continue
        si = sample_idx[row["sample"]]
        ai = asv_idx[row["sequence"]]
        total = sample_totals[row["sample"]]
        prop = round_float(cnt / total if total > 0 else 0.0, 6)
        data.append([si, ai, cnt, prop])

    result = {
        "data": data,
        "samples": samples,
        "asvs": asv_ids,
    }

    log_info(f"counts.json.gz: {len(data)} non-zero entries, "
             f"{len(samples)} samples, {len(asv_ids)} ASVs")
    return result


# ---------------------------------------------------------------------------
# Build network.json
# ---------------------------------------------------------------------------

def build_network(network_df, seq_to_id):
    """Build the network edge list.

    Returns:
        dict with 'edges' key
    """
    log_info("Building network.json")

    # Build ASV index lookup
    asv_ids_sorted = sorted(seq_to_id.values())
    asv_id_to_idx = {aid: i for i, aid in enumerate(asv_ids_sorted)}

    # Also need sequence-to-ASV-index
    seq_to_idx = {}
    for seq, aid in seq_to_id.items():
        if aid in asv_id_to_idx:
            seq_to_idx[seq] = asv_id_to_idx[aid]

    edges = []
    if network_df is not None and len(network_df) > 0:
        for _, row in network_df.iterrows():
            n1 = row["node1"]
            n2 = row["node2"]
            corr = round_float(row["correlation"])

            idx1 = seq_to_idx.get(n1)
            idx2 = seq_to_idx.get(n2)
            if idx1 is not None and idx2 is not None:
                edges.append([idx1, idx2, corr])

    log_info(f"network.json: {len(edges)} edges")
    return {"edges": edges}


# ---------------------------------------------------------------------------
# Build taxonomy.json
# ---------------------------------------------------------------------------

def build_taxonomy(taxonomy_dict, seq_to_id):
    """Build per-database taxonomy assignments.

    Returns:
        dict keyed by database name
    """
    log_info("Building taxonomy.json")

    result = {}
    for db_name, tax_data in taxonomy_dict.items():
        if "tax" in tax_data and isinstance(tax_data["tax"], pd.DataFrame):
            tax_df = tax_data["tax"]
        elif isinstance(tax_data, pd.DataFrame):
            tax_df = tax_data
        else:
            log_warn(f"Skipping taxonomy database '{db_name}': unexpected format")
            continue

        levels = list(tax_df.columns)
        assignments = {}
        for seq in tax_df.index:
            if seq in seq_to_id:
                asv_id = seq_to_id[seq]
                vals = [str(v) if pd.notna(v) else "" for v in tax_df.loc[seq]]
                assignments[asv_id] = vals

        result[db_name] = {
            "levels": levels,
            "assignments": assignments,
        }
        log_info(f"  {db_name}: {len(levels)} levels, "
                 f"{len(assignments)} ASVs assigned")

    return result


# ---------------------------------------------------------------------------
# Build renorm_stats.json
# ---------------------------------------------------------------------------

def build_renorm_stats(renorm_data):
    """Build group-level summary statistics.

    Accepts either:
      - renorm_merged DataFrame (with group, count columns)
      - renorm_by_group dict of DataFrames

    Returns:
        dict keyed by group name
    """
    log_info("Building renorm_stats.json")

    result = {}

    if isinstance(renorm_data, dict):
        # renorm_by_group format
        for grp, df in renorm_data.items():
            result[grp] = {
                "n_asvs": int(df["sequence"].nunique()) if "sequence" in df.columns else 0,
                "n_samples": int(df["sample"].nunique()) if "sample" in df.columns else 0,
                "n_reads": int(df["count"].sum()) if "count" in df.columns else 0,
            }
    elif isinstance(renorm_data, pd.DataFrame):
        # renorm_merged format
        if "group" in renorm_data.columns:
            for grp, sub in renorm_data.groupby("group"):
                result[grp] = {
                    "n_asvs": int(sub["sequence"].nunique()),
                    "n_samples": int(sub["sample"].nunique()),
                    "n_reads": int(sub["count"].sum()),
                }
        else:
            log_warn("renorm data has no 'group' column; stats will be empty")
    else:
        log_warn(f"Unexpected renorm data type: {type(renorm_data)}")

    for grp, stats in result.items():
        log_info(f"  {grp}: {stats['n_asvs']} ASVs, "
                 f"{stats['n_reads']:,} reads")

    return result


# ---------------------------------------------------------------------------
# Load taxonomy files from a directory
# ---------------------------------------------------------------------------

def load_taxonomy_dir(taxonomy_dir):
    """Load all *_taxonomy.pkl files from a directory.

    Returns:
        dict of {db_name: tax_data}
    """
    taxonomy_dict = {}

    # Search provided directory
    patterns = [
        os.path.join(taxonomy_dir, "*_taxonomy.pkl"),
        os.path.join(taxonomy_dir, "*_taxonomy.pickle"),
    ]
    tax_files = []
    for pattern in patterns:
        tax_files.extend(glob.glob(pattern))

    # Fallback: search current working directory
    if not tax_files:
        log_warn(f"No taxonomy files in {taxonomy_dir}, searching current directory")
        tax_files = glob.glob("*_taxonomy.pkl")

    if not tax_files:
        log_warn("No taxonomy files found anywhere")
        return taxonomy_dict

    for tax_file in sorted(tax_files):
        # Extract database name from filename: <db_name>_taxonomy.pkl
        basename = os.path.basename(tax_file)
        db_name = basename.replace("_taxonomy.pkl", "").replace("_taxonomy.pickle", "")
        log_info(f"Loading taxonomy: {tax_file} (db={db_name})")
        try:
            tax_data = load_pickle(tax_file)
            taxonomy_dict[db_name] = tax_data
        except Exception as e:
            log_warn(f"Failed to load {tax_file}: {e}")

    return taxonomy_dict


# ---------------------------------------------------------------------------
# Copy Svelte app
# ---------------------------------------------------------------------------

def copy_svelte_app():
    """Copy the Svelte app's index.html to the output directory if it exists."""
    # Look for the built app relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "..", "..", "viz", "dist", "index.html"),
        os.path.join(script_dir, "..", "..", "..", "viz", "dist", "index.html"),
        os.path.join(script_dir, "..", "..", "viz", "public", "index.html"),
    ]

    for candidate in candidates:
        candidate = os.path.normpath(candidate)
        if os.path.isfile(candidate):
            import shutil
            dest = os.path.join(os.getcwd(), "index.html")
            shutil.copy2(candidate, dest)
            log_info(f"Copied Svelte app from {candidate}")
            return

    log_warn("Svelte app dist/index.html not found; skipping copy")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 8:
        log_error(
            "Usage: build_viz.py <seqtab.pkl> <renorm.pkl> <taxonomy_dir> "
            "<metadata.pkl_or_NONE> <sample_tsne.pkl> <seq_tsne.pkl> <network.pkl>"
        )
        sys.exit(1)

    seqtab_path = sys.argv[1]
    renorm_path = sys.argv[2]
    taxonomy_dir = sys.argv[3]
    metadata_path = sys.argv[4]
    sample_tsne_path = sys.argv[5]
    seq_tsne_path = sys.argv[6]
    network_path = sys.argv[7]

    # --- Load data ---
    seqtab = load_pickle(seqtab_path)
    if not isinstance(seqtab, pd.DataFrame):
        log_error("seqtab.pkl must contain a pandas DataFrame")
        sys.exit(1)

    # Normalize column names
    for col in ("sample", "sequence", "count"):
        if col not in seqtab.columns:
            log_error(f"seqtab missing required column: {col}")
            sys.exit(1)

    renorm_data = load_pickle(renorm_path)

    # Extract the renorm_table_list (sequence-to-group mapping)
    renorm_table_list = None
    if isinstance(renorm_data, pd.DataFrame) and "group" in renorm_data.columns:
        # renorm_merged format: extract unique sequence->group mapping
        renorm_table_list = renorm_data[["sequence", "group"]].drop_duplicates()
    elif isinstance(renorm_data, dict):
        # Could be renorm_by_group; build the mapping from it
        rows = []
        for grp, df in renorm_data.items():
            if "sequence" in df.columns:
                for seq in df["sequence"].unique():
                    rows.append({"sequence": seq, "group": grp})
        if rows:
            renorm_table_list = pd.DataFrame(rows).drop_duplicates()

    # Load taxonomy
    taxonomy_dict = load_taxonomy_dir(taxonomy_dir)

    # Load metadata (optional)
    metadata = None
    if metadata_path and metadata_path.upper() != "NONE":
        if os.path.isfile(metadata_path):
            try:
                metadata = load_pickle(metadata_path)
                log_info(f"Metadata loaded: {type(metadata)}")
            except Exception as e:
                log_warn(f"Failed to load metadata: {e}")
        else:
            log_warn(f"Metadata file not found: {metadata_path}")

    sample_tsne = load_pickle(sample_tsne_path)
    seq_tsne = load_pickle(seq_tsne_path)

    # Load network (may be empty)
    network_df = None
    if os.path.isfile(network_path):
        network_df = load_pickle(network_path)
        if isinstance(network_df, pd.DataFrame):
            log_info(f"Network: {len(network_df)} edges")
        else:
            log_warn(f"Network file has unexpected type: {type(network_df)}")
            network_df = None

    # --- Build outputs ---

    # samples.json
    samples = build_samples(seqtab, sample_tsne, metadata)
    write_json(samples, "samples.json")

    # asvs.json.gz
    asvs, seq_to_id = build_asvs(seqtab, seq_tsne, renorm_table_list, taxonomy_dict)
    write_json(asvs, "asvs.json.gz", compress=True)

    # counts.json.gz
    counts = build_counts(seqtab, seq_to_id)
    write_json(counts, "counts.json.gz", compress=True)

    # network.json
    network = build_network(network_df, seq_to_id)
    write_json(network, "network.json")

    # taxonomy.json
    taxonomy = build_taxonomy(taxonomy_dict, seq_to_id)
    write_json(taxonomy, "taxonomy.json")

    # renorm_stats.json
    renorm_stats = build_renorm_stats(renorm_data)
    write_json(renorm_stats, "renorm_stats.json")

    # --- Copy Svelte app ---
    copy_svelte_app()

    log_info("build_viz.py complete")


if __name__ == "__main__":
    main()

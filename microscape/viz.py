"""
microscape.viz -- Export pipeline results to JSON for the Svelte frontend.

Converts DataFrames and analysis outputs into compact JSON files suitable
for client-side interactive visualization.
"""

import os
import json
import gzip
import glob as globmod
import pickle
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _round_float(val, decimals=4):
    """Round a float, returning 0.0 for NaN/Inf/None."""
    if val is None or (isinstance(val, float) and (np.isnan(val) or np.isinf(val))):
        return 0.0
    return round(float(val), decimals)


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


def _write_json(obj, path, compress=False):
    """Write an object as JSON, optionally gzip-compressed."""
    json_str = json.dumps(obj, separators=(",", ":"), allow_nan=False)

    if compress:
        gz_path = path if path.endswith(".gz") else path + ".gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
            fh.write(json_str)
    else:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(json_str)


def _load_pickle(path):
    """Load a pickle file."""
    with open(path, "rb") as fh:
        return pickle.load(fh)


# ---------------------------------------------------------------------------
# Builders (ported from build_viz.py)
# ---------------------------------------------------------------------------

def _build_samples(seqtab, sample_tsne, metadata):
    """Build the samples array with t-SNE coordinates and metadata."""
    sample_stats = seqtab.groupby("sample").agg(
        total_reads=("count", "sum"),
        n_asvs=("sequence", "nunique"),
    ).reset_index()

    tsne_lookup = {}
    for _, row in sample_tsne.iterrows():
        tsne_lookup[row["label"]] = (
            _round_float(row["tSNE1"]),
            _round_float(row["tSNE2"]),
        )

    meta_lookup = {}
    meta_fields = []
    if metadata is not None and isinstance(metadata, pd.DataFrame):
        id_col = None
        for candidate in ("sample", "Sample", "sample_id", "SampleID", "id"):
            if candidate in metadata.columns:
                id_col = candidate
                break
        if id_col is None and metadata.index.name:
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
        if sid in meta_lookup:
            rec.update(meta_lookup[sid])
        records.append(rec)

    return records


def _build_asvs(seqtab, seq_tsne, renorm_table_list, taxonomy_dict):
    """Build the ASVs array with t-SNE coordinates, group, and taxonomy."""
    sequences = sorted(seqtab["sequence"].unique())
    seq_to_id = {seq: f"ASV_{i + 1:04d}" for i, seq in enumerate(sequences)}

    asv_stats = seqtab.groupby("sequence").agg(
        total_reads=("count", "sum"),
        n_samples=("sample", "nunique"),
    ).reset_index()

    tsne_lookup = {}
    for _, row in seq_tsne.iterrows():
        tsne_lookup[row["label"]] = (
            _round_float(row["tSNE1"]),
            _round_float(row["tSNE2"]),
        )

    group_lookup = {}
    if renorm_table_list is not None:
        for _, row in renorm_table_list.iterrows():
            group_lookup[row["sequence"]] = row["group"]

    tax_string_lookup = {}
    if taxonomy_dict:
        first_db = next(iter(taxonomy_dict))
        tax_data = taxonomy_dict[first_db]
        if "tax" in tax_data and isinstance(tax_data["tax"], pd.DataFrame):
            tax_df = tax_data["tax"]
            for seq in tax_df.index:
                vals = [str(v) if pd.notna(v) else "" for v in tax_df.loc[seq]]
                tax_string_lookup[seq] = ";".join(vals)

    records = []
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

    return records, seq_to_id


def _build_counts(seqtab, seq_to_id):
    """Build the sparse count matrix."""
    samples = sorted(seqtab["sample"].unique())
    asvs = sorted(seq_to_id.keys(), key=lambda s: seq_to_id[s])
    asv_ids = [seq_to_id[s] for s in asvs]

    sample_idx = {s: i for i, s in enumerate(samples)}
    asv_idx = {s: i for i, s in enumerate(asvs)}
    sample_totals = seqtab.groupby("sample")["count"].sum()

    data = []
    for _, row in seqtab.iterrows():
        cnt = int(row["count"])
        if cnt == 0:
            continue
        si = sample_idx[row["sample"]]
        ai = asv_idx[row["sequence"]]
        total = sample_totals[row["sample"]]
        prop = _round_float(cnt / total if total > 0 else 0.0, 6)
        data.append([si, ai, cnt, prop])

    return {
        "data": data,
        "samples": samples,
        "asvs": asv_ids,
    }


def _build_network(network_df, seq_to_id):
    """Build the network edge list."""
    asv_ids_sorted = sorted(seq_to_id.values())
    asv_id_to_idx = {aid: i for i, aid in enumerate(asv_ids_sorted)}

    seq_to_idx = {}
    for seq, aid in seq_to_id.items():
        if aid in asv_id_to_idx:
            seq_to_idx[seq] = asv_id_to_idx[aid]

    edges = []
    if network_df is not None and len(network_df) > 0:
        for _, row in network_df.iterrows():
            n1 = row["node1"]
            n2 = row["node2"]
            corr = _round_float(row["correlation"])
            idx1 = seq_to_idx.get(n1)
            idx2 = seq_to_idx.get(n2)
            if idx1 is not None and idx2 is not None:
                edges.append([idx1, idx2, corr])

    return {"edges": edges}


def _build_taxonomy(taxonomy_dict, seq_to_id):
    """Build per-database taxonomy assignments."""
    result = {}
    for db_name, tax_data in taxonomy_dict.items():
        if "tax" in tax_data and isinstance(tax_data["tax"], pd.DataFrame):
            tax_df = tax_data["tax"]
        elif isinstance(tax_data, pd.DataFrame):
            tax_df = tax_data
        else:
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

    return result


def _build_renorm_stats(renorm_data):
    """Build group-level summary statistics."""
    result = {}

    if isinstance(renorm_data, dict):
        for grp, df in renorm_data.items():
            if not isinstance(df, pd.DataFrame):
                continue
            result[grp] = {
                "n_asvs": int(df["sequence"].nunique()) if "sequence" in df.columns else 0,
                "n_samples": int(df["sample"].nunique()) if "sample" in df.columns else 0,
                "n_reads": int(df["count"].sum()) if "count" in df.columns else 0,
            }
    elif isinstance(renorm_data, pd.DataFrame):
        if "group" in renorm_data.columns:
            for grp, sub in renorm_data.groupby("group"):
                result[grp] = {
                    "n_asvs": int(sub["sequence"].nunique()),
                    "n_samples": int(sub["sample"].nunique()),
                    "n_reads": int(sub["count"].sum()),
                }

    return result


def _load_taxonomy_dir(taxonomy_dir):
    """Load all ``*_taxonomy.pkl`` files from a directory."""
    taxonomy_dict = {}
    patterns = [
        os.path.join(taxonomy_dir, "*_taxonomy.pkl"),
        os.path.join(taxonomy_dir, "*_taxonomy.pickle"),
    ]
    tax_files = []
    for pattern in patterns:
        tax_files.extend(globmod.glob(pattern))

    for tax_file in sorted(tax_files):
        basename = os.path.basename(tax_file)
        db_name = basename.replace("_taxonomy.pkl", "").replace("_taxonomy.pickle", "")
        try:
            tax_data = _load_pickle(tax_file)
            taxonomy_dict[db_name] = tax_data
        except Exception:
            pass

    return taxonomy_dict


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def export_viz(results_dir, output_dir):
    """Convert pipeline outputs to JSON files for the Svelte frontend.

    Parameters
    ----------
    results_dir : str
        Directory containing pipeline pickle outputs. Expected files:

        * ``seqtab_final.pkl`` -- filtered long-format count table
        * ``renorm_merged.pkl`` or ``renorm_by_group.pkl`` -- renormalized data
        * ``sample_bray_tsne.pkl`` -- sample t-SNE coordinates
        * ``seq_bray_tsne.pkl`` -- ASV t-SNE coordinates
        * ``sparcc_correlations.pkl`` -- network edge list
        * ``metadata.pkl`` (optional) -- sample metadata
        * ``*_taxonomy.pkl`` (optional) -- taxonomy assignments

    output_dir : str
        Directory where JSON files will be written. Created if it does not
        exist.

    Notes
    -----
    Produces the following files in *output_dir*:

    * ``samples.json`` -- sample metadata with t-SNE coordinates
    * ``asvs.json.gz`` -- ASV info with t-SNE coordinates and taxonomy
    * ``counts.json.gz`` -- sparse count matrix
    * ``network.json`` -- correlation edge list
    * ``taxonomy.json`` -- per-database taxonomy assignments
    * ``renorm_stats.json`` -- group-level summary statistics
    """
    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Load inputs
    # ------------------------------------------------------------------
    seqtab_path = os.path.join(results_dir, "seqtab_final.pkl")
    seqtab = _load_pickle(seqtab_path)
    if not isinstance(seqtab, pd.DataFrame):
        raise TypeError("seqtab_final.pkl must contain a pandas DataFrame")
    for col in ("sample", "sequence", "count"):
        if col not in seqtab.columns:
            raise ValueError(f"seqtab missing required column: {col}")

    # Renormalized data
    renorm_data = None
    renorm_table_list = None
    for name in ("renorm_merged.pkl", "renorm_by_group.pkl"):
        rpath = os.path.join(results_dir, name)
        if os.path.isfile(rpath):
            renorm_data = _load_pickle(rpath)
            break

    if isinstance(renorm_data, pd.DataFrame) and "group" in renorm_data.columns:
        renorm_table_list = renorm_data[["sequence", "group"]].drop_duplicates()
    elif isinstance(renorm_data, dict):
        rows = []
        for grp, sub_df in renorm_data.items():
            if isinstance(sub_df, pd.DataFrame) and "sequence" in sub_df.columns:
                for seq in sub_df["sequence"].unique():
                    rows.append({"sequence": seq, "group": grp})
        if rows:
            renorm_table_list = pd.DataFrame(rows).drop_duplicates()

    # Taxonomy
    taxonomy_dict = _load_taxonomy_dir(results_dir)

    # Metadata (optional)
    metadata = None
    meta_path = os.path.join(results_dir, "metadata.pkl")
    if os.path.isfile(meta_path):
        try:
            metadata = _load_pickle(meta_path)
        except Exception:
            metadata = None

    # t-SNE coordinates
    sample_tsne = _load_pickle(os.path.join(results_dir, "sample_bray_tsne.pkl"))
    seq_tsne = _load_pickle(os.path.join(results_dir, "seq_bray_tsne.pkl"))

    # Network (optional)
    network_df = None
    net_path = os.path.join(results_dir, "sparcc_correlations.pkl")
    if os.path.isfile(net_path):
        network_df = _load_pickle(net_path)
        if not isinstance(network_df, pd.DataFrame):
            network_df = None

    # ------------------------------------------------------------------
    # Build and write outputs
    # ------------------------------------------------------------------

    # samples.json
    samples = _build_samples(seqtab, sample_tsne, metadata)
    _write_json(samples, os.path.join(output_dir, "samples.json"))

    # asvs.json.gz
    asvs, seq_to_id = _build_asvs(seqtab, seq_tsne, renorm_table_list, taxonomy_dict)
    _write_json(asvs, os.path.join(output_dir, "asvs.json.gz"), compress=True)

    # counts.json.gz
    counts = _build_counts(seqtab, seq_to_id)
    _write_json(counts, os.path.join(output_dir, "counts.json.gz"), compress=True)

    # network.json
    network = _build_network(network_df, seq_to_id)
    _write_json(network, os.path.join(output_dir, "network.json"))

    # taxonomy.json
    taxonomy = _build_taxonomy(taxonomy_dict, seq_to_id)
    _write_json(taxonomy, os.path.join(output_dir, "taxonomy.json"))

    # renorm_stats.json
    renorm_stats = _build_renorm_stats(renorm_data)
    _write_json(renorm_stats, os.path.join(output_dir, "renorm_stats.json"))

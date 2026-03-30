"""
microscape.metadata -- Load sample metadata from MIMARKS-compliant TSV/CSV.

Reads a metadata file, auto-detects MIMARKS fields, and optionally matches
sample identifiers against a sequence table.
"""

import os
import pandas as pd


# Standard MIMARKS fields expected in amplicon metadata
MIMARKS_FIELDS = [
    "sample_name",
    "collection_date",
    "geo_loc_name",
    "lat_lon",
    "depth",
    "env_broad_scale",
    "env_local_scale",
    "env_medium",
]


def load_metadata(path, seqtab=None, id_col="sample_name"):
    """Load sample metadata from a TSV or CSV file.

    Parameters
    ----------
    path : str
        Path to the metadata file. Files with ``.tsv`` or ``.txt``
        extensions are read as tab-separated; everything else as
        comma-separated.
    seqtab : pd.DataFrame or None
        Optional long-format sequence table (columns: sample, sequence,
        count). When provided, the returned DataFrame includes additional
        ``_match_status`` information so users can see which samples exist
        only in the pipeline output, only in the metadata, or in both.
    id_col : str
        Name of the column containing sample identifiers. If this column
        is not found, the first column of the file is used instead.

    Returns
    -------
    pd.DataFrame
        Metadata indexed by sample identifier. The DataFrame has an
        attribute ``match_stats`` (a dict) summarising the overlap between
        the metadata and the sequence table when *seqtab* is provided.
    """
    # ------------------------------------------------------------------
    # Read file
    # ------------------------------------------------------------------
    ext = os.path.splitext(path)[1].lower()
    if ext in (".tsv", ".txt"):
        meta = pd.read_csv(path, sep="\t")
    else:
        meta = pd.read_csv(path)

    # ------------------------------------------------------------------
    # Auto-detect MIMARKS fields
    # ------------------------------------------------------------------
    found_mimarks = [f for f in MIMARKS_FIELDS if f in meta.columns]

    # ------------------------------------------------------------------
    # Set index to the sample-ID column
    # ------------------------------------------------------------------
    if id_col in meta.columns:
        meta = meta.set_index(id_col)
    else:
        meta = meta.set_index(meta.columns[0])

    # ------------------------------------------------------------------
    # Match against pipeline samples when a seqtab is provided
    # ------------------------------------------------------------------
    match_stats = {
        "mimarks_fields": found_mimarks,
        "n_metadata": len(meta),
    }

    if seqtab is not None:
        if isinstance(seqtab, pd.DataFrame) and "sample" in seqtab.columns:
            pipeline_samples = set(seqtab["sample"].unique())
        else:
            pipeline_samples = set()

        matched = pipeline_samples & set(meta.index)
        in_pipe_only = pipeline_samples - set(meta.index)
        in_meta_only = set(meta.index) - pipeline_samples

        match_stats.update({
            "matched": len(matched),
            "pipeline_only": len(in_pipe_only),
            "metadata_only": len(in_meta_only),
            "total_pipeline": len(pipeline_samples),
        })

    # Attach match_stats as a DataFrame attribute so callers can inspect it
    meta.attrs["match_stats"] = match_stats

    return meta

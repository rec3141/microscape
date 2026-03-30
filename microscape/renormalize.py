"""
microscape.renormalize -- Group ASVs by taxonomy and normalize within groups.

Assigns each ASV to a biological group (prokaryote, eukaryote, chloroplast,
mitochondria, unknown) based on its taxonomic classification, then computes
within-group proportions per sample.
"""

import pandas as pd


def _find_col(df, candidates):
    """Return the first matching column name (case-insensitive)."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def _assign_group(row, col_domain, col_order, col_family):
    """Classify a sequence into a biological group based on taxonomy."""
    domain = str(row.get(col_domain, "")).strip() if col_domain else ""
    order = str(row.get(col_order, "")).strip() if col_order else ""
    family = str(row.get(col_family, "")).strip() if col_family else ""

    # Check organellar origins first (these override domain)
    if order == "Chloroplast":
        return "chloroplast"
    if family == "Mitochondria":
        return "mitochondria"
    if domain == "Eukaryota":
        return "eukaryote"
    if domain in ("Bacteria", "Archaea"):
        return "prokaryote"
    return "unknown"


def renormalize(df, taxa):
    """Group ASVs by taxonomy and compute within-group proportions.

    Parameters
    ----------
    df : pd.DataFrame
        Long-format count table with columns ``sample``, ``sequence``,
        ``count``.
    taxa : pd.DataFrame
        Taxonomy table indexed by sequence. Expected to contain columns
        for at least one of Kingdom/Domain, Order, and Family (case-
        insensitive matching is applied).

    Returns
    -------
    dict
        A dictionary with the following keys:

        * One key per biological group (e.g. ``"prokaryote"``,
          ``"eukaryote"``, ``"chloroplast"``, ``"mitochondria"``,
          ``"unknown"``), each mapping to a long-format pd.DataFrame
          with an added ``proportion`` column representing within-group
          relative abundance.
        * ``"group_assignments"`` -- a pd.DataFrame mapping each sequence
          to its assigned group.
        * ``"stats"`` -- a pd.DataFrame with per-group summary statistics.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame with columns: sample, sequence, count")
    for col in ("sample", "sequence", "count"):
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    if not isinstance(taxa, pd.DataFrame):
        raise TypeError("taxa must be a pandas DataFrame indexed by sequence")

    # ------------------------------------------------------------------
    # Identify taxonomy columns (flexible naming)
    # ------------------------------------------------------------------
    col_domain = _find_col(taxa, ["Kingdom", "Domain", "Level1"])
    col_order = _find_col(taxa, ["Order", "Level4"])
    col_family = _find_col(taxa, ["Family", "Level5"])

    # ------------------------------------------------------------------
    # Build sequence -> group mapping
    # ------------------------------------------------------------------
    seq_groups = {}
    for seq in taxa.index:
        row = taxa.loc[seq]
        seq_groups[seq] = _assign_group(row, col_domain, col_order, col_family)

    # Sequences in count table but missing from taxonomy
    all_seqs = df["sequence"].unique()
    for seq in all_seqs:
        if seq not in seq_groups:
            seq_groups[seq] = "unknown"

    # ------------------------------------------------------------------
    # Add group column and compute within-group proportions
    # ------------------------------------------------------------------
    seqtab = df.copy()
    seqtab["group"] = seqtab["sequence"].map(seq_groups)

    group_sample_totals = seqtab.groupby(["group", "sample"])["count"].transform("sum")
    seqtab["proportion"] = seqtab["count"] / group_sample_totals
    seqtab["proportion"] = seqtab["proportion"].fillna(0.0)

    # ------------------------------------------------------------------
    # Split by group
    # ------------------------------------------------------------------
    by_group = {}
    for grp, sub in seqtab.groupby("group"):
        by_group[grp] = sub.reset_index(drop=True)

    # ------------------------------------------------------------------
    # Build summary statistics
    # ------------------------------------------------------------------
    stats_rows = []
    for grp in sorted(by_group.keys()):
        sub = by_group[grp]
        stats_rows.append({
            "group": grp,
            "n_asvs": sub["sequence"].nunique(),
            "n_samples": sub["sample"].nunique(),
            "total_reads": int(sub["count"].sum()),
            "mean_reads_per_sample": round(
                sub.groupby("sample")["count"].sum().mean(), 1
            ),
            "median_reads_per_sample": round(
                sub.groupby("sample")["count"].sum().median(), 1
            ),
        })
    stats_df = pd.DataFrame(stats_rows)

    # ------------------------------------------------------------------
    # Sequence-to-group mapping table
    # ------------------------------------------------------------------
    group_assignments = pd.DataFrame({
        "sequence": list(seq_groups.keys()),
        "group": list(seq_groups.values()),
    })

    # ------------------------------------------------------------------
    # Build result
    # ------------------------------------------------------------------
    result = dict(by_group)
    result["group_assignments"] = group_assignments
    result["stats"] = stats_df
    return result

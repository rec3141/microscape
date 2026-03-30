"""
microscape.network -- SparCC-style co-occurrence network via CLR correlations.

Uses a centered log-ratio (CLR) transformation followed by Pearson correlation
to approximate SparCC, then returns an edge list filtered by prevalence and
correlation thresholds.
"""

import numpy as np
import pandas as pd


def _clr_transform(count_matrix, pseudocount=0.5):
    """Centered log-ratio transformation.

    Parameters
    ----------
    count_matrix : np.ndarray
        Matrix of shape (samples, ASVs).
    pseudocount : float
        Value added before log transformation to handle zeros.

    Returns
    -------
    np.ndarray
        CLR-transformed matrix of the same shape.
    """
    mat = count_matrix + pseudocount
    log_mat = np.log(mat)
    geo_mean_log = log_mat.mean(axis=1, keepdims=True)
    return log_mat - geo_mean_log


def sparcc_network(df, min_prevalence=0.1, min_correlation=0.1):
    """Compute a CLR-based co-occurrence network (SparCC approximation).

    Parameters
    ----------
    df : pd.DataFrame
        Long-format count table with columns ``sample``, ``sequence``,
        ``count``.
    min_prevalence : float
        Minimum fraction of samples in which an ASV must be present to be
        included in the network (0 to 1).
    min_correlation : float
        Minimum absolute Pearson correlation (on CLR-transformed data) for
        an edge to be retained in the output.

    Returns
    -------
    pd.DataFrame
        Edge list with columns ``node1``, ``node2``, ``correlation``,
        ``weight`` (absolute correlation), and ``color`` (``"positive"`` or
        ``"negative"``). Sorted by descending weight.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame with columns: sample, sequence, count")

    # Handle flexible column names
    required_cols = {"sample", "sequence", "count"}
    if not required_cols.issubset(set(df.columns)):
        col_map = {}
        for col in df.columns:
            cl = col.lower()
            if cl in ("sample", "sample_id", "sampleid"):
                col_map[col] = "sample"
            elif cl in ("sequence", "seq", "asv"):
                col_map[col] = "sequence"
            elif cl in ("count", "abundance", "reads"):
                col_map[col] = "count"
        if set(col_map.values()) >= required_cols:
            df = df.rename(columns=col_map)
        else:
            raise ValueError(
                f"Expected columns {required_cols}, got {set(df.columns)}"
            )

    # ------------------------------------------------------------------
    # Pivot to wide count matrix (samples x sequences)
    # ------------------------------------------------------------------
    wide = df.pivot_table(
        index="sample", columns="sequence", values="count", fill_value=0
    )
    n_samples = wide.shape[0]

    # ------------------------------------------------------------------
    # Filter ASVs by prevalence
    # ------------------------------------------------------------------
    prevalence = (wide > 0).sum(axis=0) / n_samples
    keep_mask = prevalence >= min_prevalence
    n_kept = int(keep_mask.sum())

    if n_kept < 2:
        return pd.DataFrame(
            columns=["node1", "node2", "correlation", "weight", "color"]
        )

    wide_filtered = wide.loc[:, keep_mask]
    asv_labels = list(wide_filtered.columns)

    # ------------------------------------------------------------------
    # CLR transform and Pearson correlation
    # ------------------------------------------------------------------
    clr_matrix = _clr_transform(wide_filtered.values)
    cor_matrix = np.corrcoef(clr_matrix.T)

    # ------------------------------------------------------------------
    # Melt to edge list (upper triangle, no self-loops)
    # ------------------------------------------------------------------
    n = len(asv_labels)
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            corr_val = cor_matrix[i, j]
            if np.isnan(corr_val):
                continue
            if abs(corr_val) >= min_correlation:
                edges.append({
                    "node1": asv_labels[i],
                    "node2": asv_labels[j],
                    "correlation": corr_val,
                })

    edge_df = pd.DataFrame(edges)

    if len(edge_df) > 0:
        edge_df["weight"] = edge_df["correlation"].abs()
        edge_df["color"] = edge_df["correlation"].apply(
            lambda x: "positive" if x > 0 else "negative"
        )
        edge_df = edge_df.sort_values("weight", ascending=False).reset_index(drop=True)
    else:
        edge_df = pd.DataFrame(
            columns=["node1", "node2", "correlation", "weight", "color"]
        )

    return edge_df

"""
microscape.ordination -- Compute ordinations (t-SNE, PCA) for samples and ASVs.

Pivots a long-format sequence table to a proportional matrix, computes pairwise
Bray-Curtis (or other) distances, and produces low-dimensional embeddings via
PCA followed by t-SNE.
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def _compute_distance_matrix(prop_matrix, metric="braycurtis"):
    """Compute a pairwise distance matrix from a proportional abundance matrix."""
    dist_vec = pdist(prop_matrix, metric=metric)
    return squareform(dist_vec)


def _run_ordination(dist_matrix, labels, method="tsne", perplexity=30.0):
    """Run dimensionality reduction on a distance matrix.

    Returns a DataFrame with columns label, dim1, dim2.
    """
    n = dist_matrix.shape[0]

    # PCA on the distance matrix as an initial reduction step
    n_pca = min(50, n - 1) if n > 1 else 1
    pca = PCA(n_components=n_pca)
    pca_coords = pca.fit_transform(dist_matrix)

    if method == "tsne" and n >= 3:
        effective_perplexity = perplexity
        if n <= perplexity:
            effective_perplexity = max(1.0, float(n - 1) / 2.0)

        tsne = TSNE(
            n_components=2,
            perplexity=effective_perplexity,
            max_iter=1000,
            random_state=42,
        )
        coords = tsne.fit_transform(pca_coords)
    else:
        # Fewer than 3 entities or PCA requested -- use first two PCA components
        dim1 = pca_coords[:, 0] if pca_coords.shape[1] >= 1 else np.zeros(n)
        dim2 = pca_coords[:, 1] if pca_coords.shape[1] >= 2 else np.zeros(n)
        coords = np.column_stack([dim1, dim2])

    return pd.DataFrame({
        "label": labels,
        "dim1": coords[:, 0],
        "dim2": coords[:, 1],
    })


def ordinate(df, method="tsne", metric="braycurtis", perplexity=30):
    """Compute ordination coordinates for samples and ASVs.

    Parameters
    ----------
    df : pd.DataFrame
        Long-format count table with columns ``sample``, ``sequence``,
        ``count``.
    method : str
        Dimensionality reduction method. Currently supports ``"tsne"``
        (PCA followed by t-SNE) and ``"pca"`` (PCA only).
    metric : str
        Distance metric passed to :func:`scipy.spatial.distance.pdist`.
        Common choices: ``"braycurtis"``, ``"euclidean"``, ``"jaccard"``.
    perplexity : float
        Perplexity parameter for t-SNE. Automatically adjusted downwards
        when the number of entities is small.

    Returns
    -------
    dict
        ``"sample_coords"``  -- pd.DataFrame with columns ``label``,
                                ``dim1``, ``dim2`` for samples.
        ``"asv_coords"``     -- pd.DataFrame with columns ``label``,
                                ``dim1``, ``dim2`` for ASVs.
        ``"distances"``      -- np.ndarray, sample-by-sample distance matrix.
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
    # Pivot to wide matrix (samples x sequences) and normalise to proportions
    # ------------------------------------------------------------------
    wide = df.pivot_table(
        index="sample", columns="sequence", values="count", fill_value=0
    )
    row_sums = wide.sum(axis=1).replace(0, 1)
    prop_matrix = wide.div(row_sums, axis=0)

    sample_labels = list(prop_matrix.index)
    seq_labels = list(prop_matrix.columns)

    # ------------------------------------------------------------------
    # Sample ordination
    # ------------------------------------------------------------------
    sample_dist = _compute_distance_matrix(prop_matrix.values, metric=metric)
    sample_coords = _run_ordination(
        sample_dist, sample_labels, method=method, perplexity=perplexity
    )

    # ------------------------------------------------------------------
    # ASV ordination (transpose)
    # ------------------------------------------------------------------
    seq_dist = _compute_distance_matrix(prop_matrix.T.values, metric=metric)
    asv_coords = _run_ordination(
        seq_dist, seq_labels, method=method, perplexity=perplexity
    )

    return {
        "sample_coords": sample_coords,
        "asv_coords": asv_coords,
        "distances": sample_dist,
    }

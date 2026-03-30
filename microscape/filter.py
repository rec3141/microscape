"""
microscape.filter -- Quality-control filtering on long-format sequence data.

Applies a cascade of filters to remove noise and low-confidence observations
from a long-format pandas DataFrame (columns: sample, sequence, count).
Operates entirely in long format -- never builds the full dense matrix, so
memory use is proportional to non-zero entries.

Filters (in order):
    1. Length      -- ASVs shorter than a threshold (can't assign taxonomy)
    2. Prevalence  -- ASVs in fewer than N samples (likely artefacts)
    3. Abundance   -- ASVs with fewer than N total reads (singletons/doubletons)
    4. Depth       -- samples with too few total reads (insufficient coverage)
"""

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def filter_seqtab(df, min_length=50, min_samples=2, min_seqs=2, min_reads=1000):
    """Apply cascading quality filters to a long-format sequence table.

    Parameters
    ----------
    df : pd.DataFrame
        Long-format DataFrame with columns ``sample``, ``sequence``, ``count``.
    min_length : int
        Minimum ASV sequence length in base-pairs. ASVs shorter than this are
        removed because they cannot be reliably assigned taxonomy.
    min_samples : int
        Minimum number of samples an ASV must appear in to be retained
        (prevalence filter). ASVs below this threshold are considered orphans.
    min_seqs : int
        Minimum total read count for an ASV across all samples (abundance
        filter). ASVs below this are considered singletons/doubletons.
    min_reads : int
        Minimum total reads per sample (depth filter). Samples with fewer
        reads are removed as insufficiently sequenced.

    Returns
    -------
    dict
        ``"filtered"``      -- pd.DataFrame of retained observations.
        ``"orphans"``       -- pd.DataFrame of ASVs removed by the prevalence
                               filter.
        ``"small_samples"`` -- pd.DataFrame of samples removed by the depth
                               filter.
        ``"stats"``         -- pd.DataFrame summarising each filter step.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Expected a pandas DataFrame with columns: sample, sequence, count")
    for col in ("sample", "sequence", "count"):
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    dt = df.copy()
    total_input = int(dt["count"].sum())
    n_input_asvs = dt["sequence"].nunique()

    # ------------------------------------------------------------------
    # 1. Remove short sequences
    # ------------------------------------------------------------------
    seq_lengths = dt.groupby("sequence").agg(
        seq_len=("sequence", lambda x: len(x.iloc[0]))
    ).reset_index()
    short_seqs = set(seq_lengths.loc[seq_lengths["seq_len"] < min_length, "sequence"])
    n_short = len(short_seqs)

    if n_short > 0:
        dt = dt[~dt["sequence"].isin(short_seqs)].copy()

    # ------------------------------------------------------------------
    # 2. Remove low-prevalence ASVs (orphans)
    # ------------------------------------------------------------------
    seq_prevalence = dt.groupby("sequence")["sample"].nunique().reset_index()
    seq_prevalence.columns = ["sequence", "n_samples"]
    orphan_seqs = set(
        seq_prevalence.loc[seq_prevalence["n_samples"] < min_samples, "sequence"]
    )
    n_orphans = len(orphan_seqs)

    dt_orphans = dt[dt["sequence"].isin(orphan_seqs)].copy()
    dt = dt[~dt["sequence"].isin(orphan_seqs)].copy()

    # ------------------------------------------------------------------
    # 3. Remove low-abundance ASVs
    # ------------------------------------------------------------------
    seq_abundance = dt.groupby("sequence")["count"].sum().reset_index()
    seq_abundance.columns = ["sequence", "total"]
    rare_seqs = set(seq_abundance.loc[seq_abundance["total"] < min_seqs, "sequence"])
    n_rare = len(rare_seqs)

    dt = dt[~dt["sequence"].isin(rare_seqs)].copy()

    # ------------------------------------------------------------------
    # 4. Remove shallow samples
    # ------------------------------------------------------------------
    sample_depth = dt.groupby("sample")["count"].sum().reset_index()
    sample_depth.columns = ["sample", "total"]
    small_samples = set(sample_depth.loc[sample_depth["total"] < min_reads, "sample"])
    n_small = len(small_samples)

    dt_small = dt[dt["sample"].isin(small_samples)].copy()
    dt = dt[~dt["sample"].isin(small_samples)].copy()

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    n_final_samples = dt["sample"].nunique()
    n_final_asvs = dt["sequence"].nunique()
    n_final_reads = int(dt["count"].sum())
    pct_retained = round(n_final_reads / max(total_input, 1) * 100, 1)

    stats = pd.DataFrame([
        {"step": "length",     "asvs_removed": n_short,   "samples_removed": None,
         "remaining_samples": None, "remaining_asvs": None, "pct_reads_retained": None},
        {"step": "prevalence", "asvs_removed": n_orphans,  "samples_removed": None,
         "remaining_samples": None, "remaining_asvs": None, "pct_reads_retained": None},
        {"step": "abundance",  "asvs_removed": n_rare,     "samples_removed": None,
         "remaining_samples": None, "remaining_asvs": None, "pct_reads_retained": None},
        {"step": "depth",      "asvs_removed": None,       "samples_removed": n_small,
         "remaining_samples": None, "remaining_asvs": None, "pct_reads_retained": None},
        {"step": "final",      "asvs_removed": None,       "samples_removed": None,
         "remaining_samples": n_final_samples, "remaining_asvs": n_final_asvs,
         "pct_reads_retained": pct_retained},
    ])

    return {
        "filtered": dt,
        "orphans": dt_orphans,
        "small_samples": dt_small,
        "stats": stats,
    }


def plot_filter_summary(df, output=None):
    """Generate diagnostic rank-abundance plots for a filtered sequence table.

    Parameters
    ----------
    df : pd.DataFrame
        Long-format DataFrame with columns ``sample``, ``sequence``, ``count``.
        Typically the ``"filtered"`` value returned by :func:`filter_seqtab`.
    output : str or None
        Path to a PDF file to save the plot. If *None*, the matplotlib Figure
        is returned instead of being saved.

    Returns
    -------
    matplotlib.figure.Figure or None
        The figure object when *output* is None; otherwise None (the figure
        is saved and closed).
    """
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    if len(df) > 0:
        reads_per_asv = df.groupby("sequence")["count"].sum().sort_values()
        reads_per_sample = df.groupby("sample")["count"].sum().sort_values()
        samps_per_asv = df.groupby("sequence")["sample"].nunique().sort_values()
        asvs_per_sample = df.groupby("sample")["sequence"].nunique().sort_values()

        ax = axes[0, 0]
        ax.plot(
            range(len(reads_per_asv)),
            np.log10(reads_per_asv.values.astype(float)),
            "o", markersize=1, color="steelblue",
        )
        ax.set_title("Reads per ASV")
        ax.set_ylabel("log10(reads)")
        ax.set_xlabel("ASV rank")

        ax = axes[0, 1]
        ax.plot(
            range(len(reads_per_sample)),
            np.log10(reads_per_sample.values.astype(float)),
            "o", markersize=1, color="steelblue",
        )
        ax.set_title("Reads per Sample")
        ax.set_ylabel("log10(reads)")
        ax.set_xlabel("Sample rank")

        ax = axes[1, 0]
        ax.plot(
            range(len(samps_per_asv)),
            np.log10(samps_per_asv.values.astype(float)),
            "o", markersize=1, color="darkgreen",
        )
        ax.set_title("Samples per ASV")
        ax.set_ylabel("log10(samples)")
        ax.set_xlabel("ASV rank")

        ax = axes[1, 1]
        ax.plot(
            range(len(asvs_per_sample)),
            np.log10(asvs_per_sample.values.astype(float)),
            "o", markersize=1, color="darkgreen",
        )
        ax.set_title("ASVs per Sample")
        ax.set_ylabel("log10(ASVs)")
        ax.set_xlabel("Sample rank")
    else:
        axes[0, 0].text(
            0.5, 0.5, "No data remaining after filtering",
            ha="center", va="center", fontsize=14,
            transform=axes[0, 0].transAxes,
        )

    fig.tight_layout()

    if output is not None:
        with PdfPages(output) as pdf:
            pdf.savefig(fig)
        plt.close("all")
        return None

    return fig

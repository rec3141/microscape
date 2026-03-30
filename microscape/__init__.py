"""
microscape -- downstream analysis tools for amplicon sequencing data.

Designed to work with papa2's output, providing filtering, metadata loading,
renormalization, phylogenetic tree construction, ordination, co-occurrence
network analysis, and visualization export.
"""

__version__ = "0.1.0"

from microscape.filter import filter_seqtab, plot_filter_summary
from microscape.metadata import load_metadata
from microscape.renormalize import renormalize
from microscape.phylogeny import build_phylogeny
from microscape.ordination import ordinate
from microscape.network import sparcc_network
from microscape.viz import export_viz

__all__ = [
    "filter_seqtab",
    "plot_filter_summary",
    "load_metadata",
    "renormalize",
    "build_phylogeny",
    "ordinate",
    "sparcc_network",
    "export_viz",
]

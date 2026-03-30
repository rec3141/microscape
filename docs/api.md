# API Reference

This page documents all public functions in microscape. Each function's
docstrings are rendered directly from the source.

---

## microscape.filter -- Sequence Table Filtering

Quality-control filtering of the ASV sequence table by length, prevalence,
abundance, and sample depth.

### filter_seqtab

::: microscape.filter.filter_seqtab
    options:
      heading_level: 3

### plot_filter_summary

::: microscape.filter.plot_filter_summary
    options:
      heading_level: 3

---

## microscape.metadata -- Metadata Loading

Load and validate MIMARKS-compliant sample metadata.

### load_metadata

::: microscape.metadata.load_metadata
    options:
      heading_level: 3

---

## microscape.renormalize -- Taxonomic Renormalization

Split ASVs into taxonomic groups and compute within-group proportional
abundances.

### renormalize

::: microscape.renormalize.renormalize
    options:
      heading_level: 3

---

## microscape.phylogeny -- Phylogenetic Tree Construction

Multiple sequence alignment (MAFFT) and neighbor-joining tree construction.

### build_phylogeny

::: microscape.phylogeny.build_phylogeny
    options:
      heading_level: 3

---

## microscape.ordination -- Ordination

Bray-Curtis dissimilarity with t-SNE or PCA dimensionality reduction.

### ordinate

::: microscape.ordination.ordinate
    options:
      heading_level: 3

---

## microscape.network -- Co-occurrence Networks

SparCC-style compositional correlation networks.

### sparcc_network

::: microscape.network.sparcc_network
    options:
      heading_level: 3

---

## microscape.viz -- Visualization Export

Export all analysis results as a JSON bundle for web-based visualization.

### export_viz

::: microscape.viz.export_viz
    options:
      heading_level: 3

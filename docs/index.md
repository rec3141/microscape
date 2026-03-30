# microscape

**Downstream analysis toolkit for amplicon sequencing data.**

microscape is a Python package and Nextflow pipeline for downstream analysis of amplicon sequencing data. It takes the output of [papa2](https://github.com/rec3141/papa2) (or any DADA2-style sequence table) and provides QC filtering, metadata integration, taxonomic renormalization, phylogenetic tree construction, ordination, co-occurrence network analysis, and visualization export.

## What does microscape do?

After DADA2 denoising produces a sequence table and taxonomy assignments, microscape handles everything that comes next:

1. **Filter** the sequence table by ASV length, prevalence, abundance, and sample depth
2. **Load metadata** from MIMARKS-compliant sample sheets
3. **Renormalize** by splitting ASVs into taxonomic groups (prokaryote, eukaryote, chloroplast, mitochondria) and computing within-group proportional abundances
4. **Build phylogeny** via MAFFT multiple sequence alignment and neighbor-joining tree construction
5. **Ordinate** samples using Bray-Curtis dissimilarity with t-SNE or PCA
6. **Compute networks** of co-occurring ASVs using SparCC correlation
7. **Export** all results as a JSON bundle for web-based visualization

## Architecture

microscape has two layers:

- **Python package** (`microscape/`) -- 8 composable functions for scripting and notebooks
- **Nextflow pipeline** (`nextflow/`) -- reproducible, scalable execution from raw reads through visualization. DADA2 steps use [papa2](https://github.com/rec3141/papa2) for a fully Python-native amplicon workflow.

## Project Layout

```
microscape/          Python package
  __init__.py          Public API (8 functions)
  filter.py            filter_seqtab, plot_filter_summary
  metadata.py          load_metadata
  renormalize.py       renormalize
  phylogeny.py         build_phylogeny
  ordination.py        ordinate
  network.py           sparcc_network
  viz.py               export_viz
nextflow/            Nextflow pipeline
  main.nf              Workflow orchestration
  nextflow.config      Parameters, profiles, resources
  modules/             Process definitions
  bin/                 Standalone R scripts (DADA2 steps)
  envs/                Conda environment specs
  primers/             Primer FASTA files
conda/               Bioconda recipe
```

## Quick Smoke Test

```python
import microscape
print(microscape.__version__)
# 0.1.0
```

## Citation

If you use microscape in published research, please cite the original DADA2 paper:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
> **DADA2: High-resolution sample inference from Illumina amplicon data.**
> *Nature Methods*, 13, 581-583. https://doi.org/10.1038/nmeth.3869

## License

BSD-3-Clause -- see [LICENSE](https://github.com/rec3141/microscape/blob/main/LICENSE).

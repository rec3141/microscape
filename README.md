# microscape

**Downstream analysis toolkit for amplicon sequencing data**

[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![Container](https://github.com/rec3141/microscape/actions/workflows/build-container.yml/badge.svg)](https://github.com/rec3141/microscape/actions/workflows/build-container.yml)

`microscape` is a Python package and Nextflow pipeline for downstream analysis of amplicon sequencing data. It takes the output of [papa2](https://github.com/rec3141/papa2) (or any DADA2-style sequence table) and provides QC filtering, metadata integration, taxonomic renormalization, phylogenetic tree construction, ordination, co-occurrence network analysis, and visualization export.

Full documentation: **https://rec3141.github.io/microscape**

---

## Key Features

- **Python-first** -- 8 composable functions covering the full post-denoising workflow
- **Nextflow pipeline** -- reproducible, scalable execution from raw reads to visualization
- **Works with papa2** -- designed as the downstream companion to [papa2](https://github.com/rec3141/papa2) for a fully Python-native amplicon workflow
- **Flexible filtering** -- length, prevalence, abundance, and depth thresholds
- **Multi-database taxonomy** -- parallel classification against SILVA, PR2, UNITE, or any custom reference
- **Ordination and networks** -- Bray-Curtis t-SNE/PCA and SparCC correlation networks
- **Container-ready** -- Docker and Apptainer images on GHCR

---

## Quick Install

### Bioconda (recommended)

```bash
conda install -c bioconda microscape
```

### From source

```bash
git clone https://github.com/rec3141/microscape.git
cd microscape
pip install -e .
```

### Container

```bash
# Docker
docker pull ghcr.io/rec3141/microscape:latest
docker run -v $(pwd):/data ghcr.io/rec3141/microscape python3 my_script.py

# Apptainer (HPC)
apptainer pull microscape.sif docker://ghcr.io/rec3141/microscape:latest
apptainer exec microscape.sif python3 my_script.py
```

---

## Usage

```python
import microscape

# 1. Filter the sequence table (length, prevalence, abundance, depth)
seqtab = microscape.filter_seqtab(
    "seqtab_nochim.pkl",
    min_seq_length=50, min_samples=2, min_seqs=3, min_reads=100,
)
microscape.plot_filter_summary(seqtab, output="filter_summary.png")

# 2. Load sample metadata (MIMARKS-compliant TSV)
meta = microscape.load_metadata("metadata.tsv", seqtab)

# 3. Renormalize by taxonomic group (prokaryote, eukaryote, chloroplast, mitochondria)
norm_tables = microscape.renormalize(seqtab, taxonomy="taxonomy_silva.pkl")

# 4. Build phylogenetic tree (MAFFT alignment + neighbor-joining)
tree = microscape.build_phylogeny(seqtab)

# 5. Ordination (Bray-Curtis t-SNE or PCA)
coords = microscape.ordinate(seqtab, method="tsne")

# 6. Co-occurrence network (SparCC)
network = microscape.sparcc_network(seqtab)

# 7. Export visualization bundle (JSON for web viewer)
microscape.export_viz(seqtab, meta, taxonomy, tree, coords, network, outdir="viz/")
```

---

## Nextflow Pipeline

The Nextflow pipeline runs the full workflow from raw reads to visualization. DADA2 denoising steps use [papa2](https://github.com/rec3141/papa2).

```bash
# Install papa2 for the denoising steps
conda install -c bioconda papa2

# Run the pipeline
nextflow run nextflow/main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -resume

# Multiple databases + phylogeny
nextflow run nextflow/main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus;pr2:/db/pr2.fasta:Domain,Supergroup,Division,Class,Order,Family,Genus,Species" \
    --run_phylogeny \
    --dada_cpus 16 \
    -resume
```

Run `nextflow run nextflow/main.nf --help` for all options.

### Pipeline Stages

#### Stage A: Preprocessing and DADA2

| Step | Process | Description |
|------|---------|-------------|
| 1 | `DEMULTIPLEX` | Optional inner-barcode demultiplexing (Mr_Demuxy) |
| 2 | `REMOVE_PRIMERS` | Primer trimming with cutadapt |
| 3 | `DADA2_FILTER_TRIM` | Per-sample quality filtering (maxEE, truncQ, PhiX removal) |
| 4 | `DADA2_LEARN_ERRORS` | Per-plate error model learning |
| 5 | `DADA2_DENOISE` | Denoising, pair merging, per-plate chimera removal |
| 6 | `MERGE_SEQTABS` | Merge per-plate tables (long-format, memory-efficient) |
| 7 | `REMOVE_CHIMERAS` | Sparse consensus chimera removal on merged data |
| 8 | `FILTER_SEQTAB` | Length, prevalence, abundance, and depth filtering |

#### Stage B: Taxonomy, Phylogeny, and Normalization

| Step | Process | Description |
|------|---------|-------------|
| 9 | `ASSIGN_TAXONOMY` | Naive Bayesian classification (one task per ref DB, parallel) |
| 10 | `BUILD_PHYLOGENY` | MAFFT alignment + NJ tree (optional, `--run_phylogeny`) |
| 11 | `RENORMALIZE` | Group ASVs by taxonomy and normalize within groups |

#### Stage C: Ordination and Networks

| Step | Process | Description |
|------|---------|-------------|
| 12 | `LOAD_METADATA` | Sample metadata integration |
| 13 | `CLUSTER_TSNE` | t-SNE ordination of samples and ASVs |
| 14 | `NETWORK_ANALYSIS` | SparCC correlation networks |
| 15 | `EXPORT_VIZ` | JSON export for web visualization |

---

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

---

## Citation

If you use microscape in published research, please cite the original DADA2 paper:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
> **DADA2: High-resolution sample inference from Illumina amplicon data.**
> *Nature Methods*, 13, 581-583. https://doi.org/10.1038/nmeth.3869

---

## License

BSD-3-Clause -- see [LICENSE](LICENSE).

# microscape

**Downstream analysis toolkit for amplicon sequencing data**

[![Docs](https://github.com/rec3141/microscape/actions/workflows/docs.yml/badge.svg)](https://rec3141.github.io/microscape/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)

`microscape` is a Python package for downstream analysis of amplicon sequencing data. It takes the output of [papa2](https://github.com/rec3141/papa2) (or any DADA2-style sequence table) and provides QC filtering, metadata integration, taxonomic renormalization, phylogenetic tree construction, ordination, co-occurrence network analysis, and visualization export.

Full documentation: **https://rec3141.github.io/microscape**

---

## Key Features

- **Python-first** -- 8 composable functions covering the full post-denoising workflow
- **Works with papa2** -- designed as the downstream companion to [papa2](https://github.com/rec3141/papa2) for a fully Python-native amplicon workflow
- **Nextflow pipeline** -- see [microscape-nf](https://github.com/rec3141/microscape-nf) for the full pipeline
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

For a complete pipeline from raw reads to visualization, see [microscape-nf](https://github.com/rec3141/microscape-nf):

```bash
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -profile conda -resume
```

---

## Related Packages

| Package | Language | Channel | Description |
|---|---|---|---|
| [papa2](https://github.com/rec3141/papa2) | Python | bioconda | DADA2 denoising |
| [microscape](https://github.com/rec3141/microscape) | Python | bioconda | Downstream analysis (this package) |
| [microscapeR](https://github.com/rec3141/microscapeR) | R | Bioconductor | R companion package |
| [microscape-nf](https://github.com/rec3141/microscape-nf) | Nextflow | GitHub | Full pipeline |

---

## Project Layout

```
microscape/          Python package
  filter.py            filter_seqtab, plot_filter_summary
  metadata.py          load_metadata
  renormalize.py       renormalize
  phylogeny.py         build_phylogeny
  ordination.py        ordinate
  network.py           sparcc_network
  viz.py               export_viz
conda/               Bioconda recipe
docs/                MkDocs documentation
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

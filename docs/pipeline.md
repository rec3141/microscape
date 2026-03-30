# Nextflow Pipeline

The microscape Nextflow pipeline runs the full amplicon analysis workflow from
raw paired-end FASTQ files to visualization output. DADA2 denoising steps use
[papa2](https://github.com/rec3141/papa2).

---

## Quick Start

```bash
# Install papa2 for the denoising steps
conda install -c bioconda papa2

# Basic run with SILVA taxonomy
nextflow run nextflow/main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva_train_set.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -resume

# Multiple databases + phylogeny
nextflow run nextflow/main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus;pr2:/db/pr2.fasta:Domain,Supergroup,Division,Class,Order,Family,Genus,Species" \
    --run_phylogeny \
    --dada_cpus 16 \
    -resume

# With persistent cache (skip completed steps across runs)
nextflow run nextflow/main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    --store_dir /scratch/microscape_cache \
    -resume
```

Run `nextflow run nextflow/main.nf --help` for all options.

---

## Pipeline Stages

### Stage A: Preprocessing and DADA2

| Step | Process | Description |
|------|---------|-------------|
| 1 | `DEMULTIPLEX` | Optional inner-barcode demultiplexing (Mr_Demuxy) |
| 2 | `REMOVE_PRIMERS` | Primer trimming with cutadapt (auto-selects by 16S/18S/ITS prefix) |
| 3 | `DADA2_FILTER_TRIM` | Per-sample quality filtering (maxEE, truncQ, PhiX removal) |
| 4 | `DADA2_LEARN_ERRORS` | Per-plate error model learning (plates share PCR history) |
| 5 | `DADA2_DENOISE` | Denoising, pair merging, per-plate chimera removal |
| 6 | `MERGE_SEQTABS` | Merge per-plate tables (long-format, memory-efficient) |
| 7 | `REMOVE_CHIMERAS` | Sparse consensus chimera removal on merged data |
| 8 | `FILTER_SEQTAB` | Length, prevalence, abundance, and depth filtering |

### Stage B: Taxonomy, Phylogeny, and Normalization

| Step | Process | Description |
|------|---------|-------------|
| 9 | `ASSIGN_TAXONOMY` | Naive Bayesian classification (one task per ref DB, parallel) |
| 10 | `BUILD_PHYLOGENY` | MAFFT alignment + NJ tree (optional, `--run_phylogeny`) |
| 11 | `RENORMALIZE` | Group ASVs by taxonomy and normalize within groups |

### Stage C: Ordination and Networks

| Step | Process | Description |
|------|---------|-------------|
| 12 | `LOAD_METADATA` | Sample metadata integration |
| 13 | `CLUSTER_TSNE` | t-SNE ordination of samples and ASVs |
| 14 | `NETWORK_ANALYSIS` | SparCC correlation networks |
| 15 | `EXPORT_VIZ` | JSON export for web visualization |

---

## Parameters

### Input/Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Directory of paired-end `*.fastq.gz` files |
| `--outdir` | `results` | Output directory |
| `--store_dir` | off | Persistent cache (skip completed steps across runs) |

### Primer Removal

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--primer_auto` | `true` | Auto-select primer file by filename prefix |
| `--primers` | auto | Override with a specific primer FASTA |
| `--primer_error_rate` | `0.12` | Cutadapt max error rate |

### DADA2

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--maxEE` | `2` | Max expected errors per read |
| `--truncQ` | `11` | Truncate at first base with quality <= Q |
| `--maxN` | `0` | Max Ns allowed |
| `--truncLen_fwd` | `0` | Truncate forward reads at position N (0 = off) |
| `--truncLen_rev` | `0` | Truncate reverse reads at position N (0 = off) |
| `--min_overlap` | `10` | Min overlap for pair merging |

### QC Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_seq_length` | `50` | Remove ASVs shorter than N bp |
| `--min_samples` | `2` | Remove ASVs in fewer than N samples |
| `--min_seqs` | `3` | Remove ASVs with fewer than N total reads |
| `--min_reads` | `100` | Remove samples with fewer than N reads |

### Taxonomy

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ref_databases` | none | Ref DBs (`"name:path:Levels;..."`) |
| `--run_phylogeny` | `false` | Build phylogenetic tree |

### Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threads` | `8` | General thread count |
| `--dada_cpus` | `8` | CPUs for DADA2 processes |
| `--dada_memory` | `16 GB` | Memory for DADA2 processes |

---

## Profiles

```bash
-profile standard   # Local execution (default)
-profile test       # Reduced resources for testing
-profile slurm      # Submit to SLURM cluster
```

---

## Data Format

The pipeline uses **long-format** tables as its canonical representation
from the merge step onward:

```
sample          sequence            count
plate1_A01      TACGGAGGATGCGA...   1523
plate1_A01      TACGGAGGATCCGA...   847
plate1_A02      TACGGAGGATGCGA...   2041
```

This avoids the memory cost of materializing a sparse dense matrix (samples x
ASVs). For large datasets (4K+ samples, 100K+ ASVs), the dense matrix can
exceed available memory while the long format uses only the non-zero entries.

---

## Dependencies

- [Nextflow](https://nextflow.io/) >= 23.04.0
- [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) (environments created automatically)
- [papa2](https://github.com/rec3141/papa2) for DADA2 denoising steps
- [cutadapt](https://cutadapt.readthedocs.io/) for primer removal

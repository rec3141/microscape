# Quickstart

## Installation

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

### Verify

```python
import microscape
print(microscape.__version__)
# 0.1.0
```

---

## Minimal Example

This example assumes you already have a sequence table and taxonomy from
[papa2](https://github.com/rec3141/papa2) or R's DADA2. Install papa2 if you
need the denoising steps:

```bash
conda install -c bioconda papa2
```

### Post-denoising workflow

```python
import microscape

# 1. Filter the sequence table
seqtab = microscape.filter_seqtab(
    "seqtab_nochim.pkl",
    min_seq_length=50,
    min_samples=2,
    min_seqs=3,
    min_reads=100,
)
microscape.plot_filter_summary(seqtab, output="filter_summary.png")

# 2. Load sample metadata
meta = microscape.load_metadata("metadata.tsv", seqtab)

# 3. Renormalize by taxonomic group
norm_tables = microscape.renormalize(seqtab, taxonomy="taxonomy_silva.pkl")

# 4. Build phylogenetic tree
tree = microscape.build_phylogeny(seqtab)

# 5. Ordination
coords = microscape.ordinate(seqtab, method="tsne")

# 6. Co-occurrence network
network = microscape.sparcc_network(seqtab)

# 7. Export visualization bundle
microscape.export_viz(seqtab, meta, norm_tables, tree, coords, network, outdir="viz/")
```

!!! note "Input formats"
    `filter_seqtab` accepts pickle (`.pkl`) or CSV files containing a
    long-format sequence table with columns `sample`, `sequence`, and `count`.
    Taxonomy files follow the same convention.

---

## Full pipeline with papa2

For a complete workflow starting from raw FASTQ files, combine papa2 for
denoising with microscape for downstream analysis:

```python
import papa2
import microscape

# --- papa2: denoising ---
fwd = ["sample1_R1.fastq.gz", "sample2_R1.fastq.gz"]
rev = ["sample1_R2.fastq.gz", "sample2_R2.fastq.gz"]

filt_fwd = ["filtered/sample1_R1.fastq.gz", "filtered/sample2_R1.fastq.gz"]
filt_rev = ["filtered/sample1_R2.fastq.gz", "filtered/sample2_R2.fastq.gz"]

papa2.filter_and_trim(fwd, filt_fwd, rev=rev, filt_rev=filt_rev,
                      trunc_len=(240, 200), max_ee=(2, 2), rm_phix=True)

errF = papa2.learn_errors(filt_fwd)
errR = papa2.learn_errors(filt_rev)

derepFs = [papa2.derep_fastq(f) for f in filt_fwd]
derepRs = [papa2.derep_fastq(f) for f in filt_rev]

dadaFs = papa2.dada(derepFs, err=errF)
dadaRs = papa2.dada(derepRs, err=errR)

mergers = [papa2.merge_pairs(dF, drF, dR, drR)
           for dF, drF, dR, drR in zip(dadaFs, derepFs, dadaRs, derepRs)]

seqtab = papa2.make_sequence_table(mergers)
seqtab_nochim = papa2.remove_bimera_denovo(seqtab)
taxa = papa2.assign_taxonomy(seqtab_nochim["seqs"], "silva_nr99_v138.1_train_set.fa.gz")

# --- microscape: downstream ---
filtered = microscape.filter_seqtab(seqtab_nochim, min_samples=2, min_seqs=3)
meta = microscape.load_metadata("metadata.tsv", filtered)
norm = microscape.renormalize(filtered, taxonomy=taxa)
tree = microscape.build_phylogeny(filtered)
coords = microscape.ordinate(filtered, method="tsne")
network = microscape.sparcc_network(filtered)
microscape.export_viz(filtered, meta, norm, tree, coords, network, outdir="viz/")
```

---

## Nextflow pipeline

For production use, the Nextflow pipeline automates the entire workflow:

```bash
nextflow run nextflow/main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -resume
```

See the [Pipeline documentation](pipeline.md) for full details.

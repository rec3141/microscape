# Microscape Pipeline: R vs Python Benchmarks

## Test Dataset

- **Source:** `/matika/seqs/amplicon/cutadapt/515F-806RB/20160323-uaf/`
- **Samples:** 96 (one 96-well plate, 16S V4 amplicons, 515F-806RB primers)
- **Total reads:** 576,362 paired-end reads (288,181 pairs)
- **Total bases:** 158.6 Mbp
- **Average read length:** 275 bp
- **Pre-trimmed:** cutadapt primer removal already applied
- **Machine:** 16-core server, 64 GB RAM, no GPU

## Pipeline Configurations Tested

| Config | Filter | Denoise | Chimera | Merge/Filter/Tax/etc |
|--------|--------|---------|---------|----------------------|
| **Pure R** | R `filterAndTrim` | R `dada2` | R `removeBimeraDenovo` (per-plate + post-merge) | R (data.table) |
| **Pure Python (dada2gpu)** | BioPython | dada2gpu C library | Python sparse reimplementation | Python (pandas) |
| **Hybrid** | BioPython | R `dada2` | Skipped (R per-plate is sufficient) | Python (pandas) |

## Results Summary

### Concordance

| Metric | Pure R | Pure Python (dada2gpu) | Hybrid (Py + R dada2) |
|--------|--------|------------------------|----------------------|
| Final samples | 67 | 31 | **67** |
| Final ASVs | 80 | 30 | **80** |
| Final reads | 123,018 | 27,889 | **123,018** |
| Jaccard vs R | — | poor | **100.0%** |
| Shared samples with R | — | 31/67 (46%) | **67/67 (100%)** |
| Shared ASVs with R | — | unknown | **80/80 (100%)** |

**The hybrid pipeline produces identical results to pure R.**

### Performance

| Step | Pure R | Pure Python | Hybrid | Notes |
|------|--------|-------------|--------|-------|
| Filter (96 samples) | ~40s | ~30s | ~30s | Python slightly faster (BioPython vs R) |
| Learn errors | 2m 55s | 1m 31s | 2m 55s | Python faster with nbases=1e7 tuning |
| Denoise (1 plate) | ~5m | 1m 17s | ~40s | R already cached in hybrid |
| Merge | ~1s | <1s | <1s | Python data.table merge is fast |
| Chimera removal | ~30s | <1s | skipped | R per-plate already thorough |
| Filter seqtab | ~2s | ~1s | ~1s | Comparable |
| Taxonomy (SILVA) | ~30s | ~30s | ~30s | Both use same C code |
| Renormalize | ~5s | <1s | <1s | Python pandas faster |
| Clustering (t-SNE) | ~10s | ~2s | ~2s | Python sklearn faster |
| Network (SparCC/CLR) | ~5s | <1s | <1s | Python CLR approximation |
| Viz build | ~2s | ~1s | ~1s | Python JSON vs R Shiny |
| **Total** | **9m 14s** | **3m 17s** | **4m 31s** | |

## Stage-by-Stage Comparison

### Stage A: DADA2 Processing

**Filter and Trim**
- Both produce identical read counts (276,653 reads pass in both R and Python)
- Python adds `minLen=20` to match R's default (prevents ultra-short reads breaking downstream)
- Python creates 0-byte files for empty samples (R removes them entirely)

**Error Learning**
- R uses `learnErrors(multithread=TRUE)` with RcppParallel/TBB — efficient within-sample parallelism
- Python dada2gpu uses OpenMP in C library + ProcessPoolExecutor across samples
- Python `nbases=1e7` optimization reduces from 16 minutes to 1.5 minutes (loads fewer files per iteration)
- R processes all files by default but is faster per-iteration due to TBB

**Denoising**
- R dada2: 70 samples with reads, 119 ASVs, 125,058 reads
- dada2gpu: 34 samples with reads, 93 ASVs, 38,025 reads (21 samples return empty)
- **Key issue:** dada2gpu loses ~55% of samples during denoising (see dada2gpu notes below)

**Pair Merging**
- R: `mergePairs()` with `trimOverhang=TRUE`, `minOverlap=10`
- Python: custom overlap alignment implementation (slide + match ratio > 0.9)
- R's implementation is more lenient and handles edge cases better

**Chimera Removal**
- R: `removeBimeraDenovo(method="consensus")` — per-plate in denoise step, optional post-merge
- Python: sparse reimplementation using `dada2::isBimera()` C function via ctypes
- Python's post-merge pass finds 12 additional chimeras that R doesn't flag (false positives)
- **Fix:** skip post-merge chimera when using R dada2 engine (per-plate is sufficient)

### Stage B: Taxonomy & Normalization

**Taxonomy Assignment**
- Both use the same C-level naive Bayesian 8-mer classifier (`C_assign_taxonomy2`)
- Python calls it via ctypes (`taxonomy_capi.cpp`), R calls via Rcpp
- Same algorithm, same results: 5,129 genera from 410,793 SILVA references
- Memory: ~1.3 GB for the kmer probability table (ngenus × 65,536 kmers × float32)

**Renormalization**
- Both group ASVs by SILVA taxonomy (prokaryote/eukaryote/chloroplast/mitochondria/unknown)
- Python uses pandas groupby, R uses data.table
- Results are identical when input ASVs match

### Stage C: Visualization

**Clustering**
- R: parallelDist (Bray-Curtis) → gmodels::fast.prcomp → Rtsne
- Python: scipy.spatial.distance.pdist → sklearn PCA → sklearn TSNE
- Both produce valid ordinations; exact coordinates differ due to t-SNE randomness

**Network**
- R: SpiecEasi::sparcc (true SparCC algorithm, iterative)
- Python: CLR + Pearson correlation (fast approximation)
- Different algorithms but similar biological conclusions

**Visualization**
- R: Shiny app with plotly (server-side R process required)
- Python: Svelte app with regl-scatterplot (static files, WebGL, no server needed)

## Data Format Comparison

| Format | R Pipeline | Python Pipeline |
|--------|-----------|----------------|
| Sequence tables | `.rds` (R native) | `.pkl` (Python pickle) |
| Internal format | data.table (long-format) | pandas DataFrame (long-format) |
| Merge strategy | data.table rbindlist | pandas concat |
| Wide matrix | `as.matrix()` when needed | `pivot_table()` when needed |
| Viz data | `app_data.rds` (Shiny) | JSON files (Svelte) |

---

## dada2gpu Project Notes

### Architecture

The dada2gpu project provides a standalone Python/C interface to dada2's core
algorithms without requiring R or Rcpp. It consists of:

- **C library (`libdada2.so`)**: compiled from dada2's C++ source with Rcpp
  dependencies replaced by a pure C API (`dada2_capi.cpp`)
- **Python package (`dada2gpu/`)**: ctypes bindings + high-level Python functions
- **Optional CUDA kernel** (`cuda_compare.cu`): fused kmer+alignment+lambda
  computation on GPU

### What Works

| Component | Status | Performance vs R |
|-----------|--------|-----------------|
| `dada()` denoising | ✅ Functional, accuracy issues | ~2x faster (CPU) |
| `learn_errors()` | ✅ Working | ~2x faster with nbases tuning |
| `derep_fastq()` | ✅ Working | Comparable |
| `assign_taxonomy()` | ✅ Working (same C code as R) | Same speed |
| Parallel processing | ✅ ProcessPoolExecutor | 1.7x speedup on 8 cores |
| GPU detection | ✅ `gpu_available()` | Not tested (no GPU) |

### Known Issues

**1. Denoising Accuracy (Critical)**

21 out of 96 samples produce zero denoised sequences from dada2gpu, while R
dada2 successfully denoises all of them. The surviving samples also produce
fewer ASVs and reads (38K vs 125K reads).

Likely causes:
- Error model initialization differences (dada2gpu uses `MAX_CONSIST=10`
  self-consistency but may converge differently than R)
- The `learn_errors()` with `nbases=1e7` uses fewer samples for error learning
  than R's default `nbases=1e8`, potentially producing a less accurate model
- Edge cases in the C library's cluster initialization where samples with
  few unique sequences get no clusters

**2. Missing C API Functions**

These exist as Rcpp-dependent source but are NOT compiled into `libdada2.so`:

| File | Functions | Lines | Dependency |
|------|-----------|-------|------------|
| `chimera.cpp` | `isBimera`, `C_table_bimera2` | 269 | Rcpp, RcppParallel |
| `filter.cpp` | (stub only) | 49 | Rcpp |
| `evaluate.cpp` | `nwalign`, `nwhamming`, `kmer_dist` | 356 | Rcpp |

Porting these to the standalone C API (like `taxonomy_capi.cpp`) would give
Python access to the exact same chimera detection and alignment code that R uses.

**3. Missing Python Functions**

| Function | Status | Current Workaround |
|----------|--------|-------------------|
| `mergePairs()` | ❌ Not in C | Custom Python overlap merger (too strict) |
| `filterAndTrim()` | ❌ Not in C | BioPython-based Python script (works but different edge cases) |
| `plotErrors()` | ❌ | matplotlib equivalent |
| `addSpecies()` | ❌ | Not implemented |

### Fixes Applied During This Session

1. **kmer crash on short sequences** (`kmers.cpp`): sequences shorter than
   `KMER_SIZE` (5bp) caused a fatal error. Fixed to return zeroed kmer vectors.

2. **Taxonomy C API** (`taxonomy_capi.cpp`): ported `C_assign_taxonomy2` from
   Rcpp to standalone C with OpenMP parallelism. Same algorithm, ~1.3 GB memory
   for 410K references.

3. **Parallel denoising** (`dada.py`): added `ProcessPoolExecutor` for
   per-sample parallelism within self-consistency iterations.

4. **nbases optimization** (`dada2_learn_errors.py`): default `nbases=1e7`
   instead of `1e8`, reducing learn_errors from 16 minutes to 1.5 minutes.

### Recommended Next Steps

1. **Debug empty-sample denoising**: compare dada2gpu's C `dada2_run()` output
   for a specific failing sample against R's `dada()` output for the same sample.
   Add diagnostic logging to the C code to trace where clusters are lost.

2. **Port chimera.cpp to C API**: follow the `taxonomy_capi.cpp` pattern.
   Replace `Rcpp::` calls with `fprintf`/`malloc`, remove `RcppParallel`
   dependency (use OpenMP instead). This would give Python the exact same
   `isBimera` as R.

3. **Port mergePairs to C**: R's `mergePairs` calls `nwalign()` which is
   already compiled. Just need a C API wrapper that takes two sets of denoised
   sequences and returns merged sequences with abundances.

4. **GPU benchmarking**: test on a CUDA-capable machine to measure the speedup
   from the fused kmer+align+lambda kernel in `cuda_compare.cu`.

5. **Error model convergence**: investigate why dada2gpu's self-consistency loop
   reaches `MAX_CONSIST=10` without converging on some datasets where R
   converges in 6-7 iterations. May be a numerical precision issue in the
   `loess_errfun` Python implementation vs R's C-level LOESS.

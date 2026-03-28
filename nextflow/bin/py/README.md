# Python implementations of microscape pipeline scripts

These scripts mirror the R implementations in `bin/` but use Python
with python-dada2 (dada2gpu) for the core denoising steps.

Select with `--lang python` when running the pipeline.

## Dependencies
- dada2gpu (local, /data/dada2_gpu)
- numpy, pandas
- cutadapt (shared with R pipeline)
- scikit-bio (for Bray-Curtis, t-SNE)
- scipy

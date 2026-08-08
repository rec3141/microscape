# Pipeline scripts

The Python engine of the amplicon workflow as it ran in the Nextflow pipeline,
preserved here when that workflow consolidated into the `illumina_amplicon`
stage of [danaSeq](https://github.com/rec3141/danaSeq).

These are command-line scripts, not package API — each takes positional
arguments and writes files. The curated equivalents are the package modules
(`microscape.filter`, `microscape.metadata`, `microscape.network`,
`microscape.ordination`, `microscape.phylogeny`, `microscape.renormalize`).

`auto_trim.py` is not duplicated here: it is `microscape.quality`, reachable as
`microscape auto-trim`, and it is the one command the pipeline ever used from
this package.

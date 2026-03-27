// Merge sequence tables and quality control.
//
// Three-step post-DADA2 processing:
//   1. MERGE_SEQTABS   — combine per-plate sequence tables (long-format output)
//   2. REMOVE_CHIMERAS — sparse consensus chimera removal (long-format I/O)
//   3. FILTER_SEQTAB   — length, prevalence, abundance filtering (long → long + wide)
//
// Data flows as long-format data.table (sample, sequence, count) through the
// entire merge → chimera → filter chain. The dense matrix is only materialized
// once at the end by FILTER_SEQTAB (seqtab_final_wide.rds) for downstream
// tools that require it.
//
// All R logic lives in bin/ scripts; these processes just wire up I/O.

process MERGE_SEQTABS {
    tag "merge-all"
    label 'process_medium'
    conda "${projectDir}/conda-envs/microscape-dada2"

    input:
    path(seqtab_files)

    output:
    path("seqtab_merged.rds"), emit: seqtab
    path("merge_stats.tsv"), emit: stats

    script:
    """
    merge_seqtabs.R
    """
}

// Sparse consensus chimera removal — reimplementation of dada2's
// removeBimeraDenovo("consensus") that operates on long-format data.
// Uses dada2::isBimera() for C-level alignment, but avoids the dense matrix.
process REMOVE_CHIMERAS {
    tag "chimera-removal"
    label 'process_high'
    conda "${projectDir}/conda-envs/microscape-dada2"

    input:
    path(seqtab_rds)

    output:
    path("seqtab_nochim.rds"), emit: seqtab
    path("chimera_stats.tsv"), emit: stats

    script:
    """
    remove_chimeras.R "${seqtab_rds}" ${task.cpus}
    """
}

// Length, prevalence, and abundance filtering of the chimera-free data.
// Produces both long-format (primary) and wide matrix (compatibility) outputs.
process FILTER_SEQTAB {
    tag "filter-qc"
    label 'process_medium'
    conda "${projectDir}/conda-envs/microscape-dada2"
    publishDir "${params.outdir}/seqtab_final", mode: 'copy'

    input:
    path(seqtab_rds)

    output:
    path("seqtab_final.rds"), emit: seqtab
    path("seqtab_final_wide.rds"), emit: seqtab_wide
    path("seqtab_orphans.rds"), emit: orphans
    path("seqtab_small.rds"), emit: small_samples
    path("filter_stats.tsv"), emit: stats
    path("sequence_summaries.pdf"), emit: plots

    script:
    """
    filter_seqtab.R \
        "${seqtab_rds}" \
        ${params.min_seq_length} ${params.min_samples} \
        ${params.min_seqs} ${params.min_reads}
    """
}

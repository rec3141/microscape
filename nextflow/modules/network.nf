// SparCC co-occurrence network analysis.
//
// Runs SparCC (Sparse Correlations for Compositional data) on the
// normalized count table to infer microbial co-occurrence networks.
// Low-prevalence ASVs are filtered out to reduce noise and computation.

process NETWORK_SPARCC {
    tag "sparcc"
    label 'process_high'
    conda "${projectDir}/conda-envs/microscape-network"
    publishDir "${params.outdir}/network", mode: 'copy'

    input:
    path(renorm_merged_rds)
    val(min_prevalence)

    output:
    path("sparcc_correlations.rds"), emit: correlations
    path("sparcc_stats.tsv"), emit: stats

    script:
    """
    network_sparcc.R "${renorm_merged_rds}" ${min_prevalence} ${task.cpus}
    """
}

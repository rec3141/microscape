// Taxonomy-based renormalization.
//
// Separates ASVs into biologically meaningful groups (prokaryotes,
// eukaryotes, chloroplasts, mitochondria, unknown) based on their
// taxonomic assignment, then normalizes counts to proportions within
// each group. This prevents organelle reads from distorting prokaryote
// relative abundance estimates.
//
// Requires a primary taxonomy assignment (typically SILVA) that provides
// Domain/Order/Family ranks for group classification.

process RENORMALIZE {
    tag "renormalize"
    label 'process_medium'
    conda "${projectDir}/conda-envs/microscape-dada2"
    publishDir "${params.outdir}/renormalized", mode: 'copy'

    input:
    path(seqtab_rds)
    tuple val(db_name), path(taxonomy_rds), path(bootstrap_rds)

    output:
    path("renorm_by_group.rds"), emit: by_group
    path("renorm_merged.rds"), emit: merged
    path("renorm_table_list.rds"), emit: table_list
    path("renorm_stats.tsv"), emit: stats

    script:
    """
    renormalize.R "${seqtab_rds}" "${taxonomy_rds}" "${db_name}"
    """
}

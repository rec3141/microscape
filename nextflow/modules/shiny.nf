// Build Shiny app data package from all pipeline outputs.
//
// Bundles sequence tables, taxonomy, t-SNE coordinates, network edges,
// and metadata into a single app_data.rds for the interactive Shiny explorer.
// Also copies app.R into the output directory so the result is self-contained.

process BUILD_SHINY {
    tag "build_shiny"
    label 'process_low'
    conda "${projectDir}/conda-envs/microscape-shiny"
    publishDir "${params.outdir}/shiny", mode: 'copy'

    input:
    path(seqtab_rds)
    path(renorm_rds)
    path(taxonomy_dir)
    path(metadata_rds)
    path(sample_tsne_rds)
    path(seq_tsne_rds)
    path(network_rds)

    output:
    path("app_data.rds"), emit: app_data
    path("app.R"),        emit: app_script

    script:
    """
    build_shiny.R \
        "${seqtab_rds}" \
        "${renorm_rds}" \
        "${taxonomy_dir}" \
        "${metadata_rds}" \
        "${sample_tsne_rds}" \
        "${seq_tsne_rds}" \
        "${network_rds}"
    """
}

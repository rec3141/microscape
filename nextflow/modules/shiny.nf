// Build visualization data package from all pipeline outputs.
//
// R mode: builds app_data.rds + app.R (Shiny)
// Python mode: builds JSON files for Svelte app (samples.json, asvs.json.gz, etc.)

process BUILD_VIZ {
    tag "build_viz"
    label 'process_low'
    conda params.lang == 'python' ? "${projectDir}/conda-envs/microscape-python" : "${projectDir}/conda-envs/microscape-r"
    publishDir "${params.outdir}/viz", mode: 'copy'

    input:
    path(seqtab)
    path(renorm)
    path(taxonomy_files)
    path(metadata)
    path(sample_tsne)
    path(seq_tsne)
    path(network)

    output:
    path("*.json"),    emit: json, optional: true
    path("*.json.gz"), emit: json_gz, optional: true
    path("app_data.rds"), emit: app_data, optional: true
    path("app.R"),        emit: app_script, optional: true

    script:
    if (params.lang == 'python')
    """
    build_viz.py \
        "${seqtab}" \
        "${renorm}" \
        "." \
        "${metadata}" \
        "${sample_tsne}" \
        "${seq_tsne}" \
        "${network}"
    """
    else
    """
    build_shiny.R \
        "${seqtab}" \
        "${renorm}" \
        "." \
        "${metadata}" \
        "${sample_tsne}" \
        "${seq_tsne}" \
        "${network}"
    """
}

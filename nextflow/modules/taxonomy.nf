// Taxonomy assignment — naive Bayesian kmer classifier.
//
// Each reference database runs as an independent task, enabling parallel
// classification across databases. Supports --lang R (dada2) or --lang python
// (custom Wang et al. 2007 implementation with multiprocessing).

def taxExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process ASSIGN_TAXONOMY {
    tag "${db_name}"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/conda-envs/microscape-python" : "${projectDir}/conda-envs/microscape-r"
    publishDir "${params.outdir}/taxonomy", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/taxonomy" : null

    input:
    tuple val(db_name), path(ref_db), val(tax_levels)
    path(seqtab)

    output:
    tuple val(db_name), path("${db_name}_taxonomy.${taxExt()}"), path("${db_name}_bootstrap.${taxExt()}"), emit: taxonomy
    path("${db_name}_taxonomy.tsv"), emit: taxonomy_tsv

    script:
    def levels_arg = tax_levels ? "\"${tax_levels}\"" : "null"
    if (params.lang == 'python')
    """
    PYTHONPATH=${params.dada2_path}:\${PYTHONPATH:-} \
    assign_taxonomy.py \
        "${seqtab}" "${ref_db}" "${db_name}" \
        ${task.cpus} ${levels_arg}
    """
    else
    """
    assign_taxonomy.R \
        "${seqtab}" "${ref_db}" "${db_name}" \
        ${task.cpus} ${levels_arg}
    """
}

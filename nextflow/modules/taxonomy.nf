// Taxonomy assignment using dada2::assignTaxonomy().
//
// Each reference database runs as an independent task, enabling parallel
// classification across databases. The naive Bayesian classifier is CPU-
// intensive but embarrassingly parallel across databases.
//
// Input is the long-format sequence table — only unique sequences are
// extracted for classification, so the count data is not needed here.

process ASSIGN_TAXONOMY {
    tag "${db_name}"
    label 'process_high'
    conda "${projectDir}/conda-envs/microscape-r"
    publishDir "${params.outdir}/taxonomy", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/taxonomy" : null

    input:
    tuple val(db_name), path(ref_db), val(tax_levels)
    path(seqtab_rds)

    output:
    tuple val(db_name), path("${db_name}_taxonomy.rds"), path("${db_name}_bootstrap.rds"), emit: taxonomy
    path("${db_name}_taxonomy.tsv"), emit: taxonomy_tsv

    script:
    def levels_arg = tax_levels ? "\"${tax_levels}\"" : "null"
    """
    assign_taxonomy.R \
        "${seqtab_rds}" "${ref_db}" "${db_name}" \
        ${task.cpus} ${levels_arg}
    """
}

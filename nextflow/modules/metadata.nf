// Load sample metadata and merge with sequence data.
//
// Reads a MIMARKS-compliant (or custom) TSV/CSV metadata file and matches
// its rows to sample IDs in the pipeline's sequence table. Produces a
// merged metadata object and matching statistics for QC review.

process LOAD_METADATA {
    tag "metadata"
    label 'process_low'
    conda "${projectDir}/conda-envs/microscape-dada2"

    input:
    path(seqtab_rds)
    path(metadata_file)
    val(sample_id_column)

    output:
    path("metadata.rds"), emit: metadata
    path("match_stats.tsv"), emit: stats

    script:
    """
    load_metadata.R "${seqtab_rds}" "${metadata_file}" "${sample_id_column}"
    """
}

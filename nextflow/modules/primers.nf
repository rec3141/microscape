// Primer removal with cutadapt.
//
// Removes primer sequences from paired-end reads. Supports auto-detection
// of primer set based on filename prefix (16S, 18S, ITS) or a user-supplied
// primer FASTA. Reads without primer matches are discarded (--discard-untrimmed).

process REMOVE_PRIMERS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/microscape-python"
    publishDir "${params.outdir}/trimmed", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/trimmed" : null

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"), path("${meta.id}_R2.trimmed.fastq.gz"), emit: reads
    path("${meta.id}_cutadapt.log"), emit: log

    script:
    // Auto-select primer file based on sample name prefix
    def primer_file = params.primers ?: (
        params.primer_auto ? (
            meta.id =~ /^16/ ? params.primers_bac :
            meta.id =~ /^18/ ? params.primers_euk :
            meta.id =~ /^ITS/ ? params.primers_its :
            params.primers_all
        ) : params.primers_all
    )
    """
    cutadapt \\
        -g file:${primer_file} \\
        -G file:${primer_file} \\
        --discard-untrimmed \\
        -j ${task.cpus} \\
        -e ${params.primer_error_rate} \\
        -o ${meta.id}_R1.trimmed.fastq.gz \\
        -p ${meta.id}_R2.trimmed.fastq.gz \\
        ${r1} ${r2} \\
        > ${meta.id}_cutadapt.log 2>&1

    echo "[INFO] ${meta.id}: primer removal complete" >&2
    """
}

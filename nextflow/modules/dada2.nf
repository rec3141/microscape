// DADA2 amplicon denoising processes.
//
// Three-step DADA2 workflow:
//   1. DADA2_FILTER_TRIM   — per-sample quality filtering
//   2. DADA2_LEARN_ERRORS  — per-plate error model learning
//   3. DADA2_DENOISE       — per-plate denoising + pair merging → sequence table
//
// Supports --lang R (default) or --lang python (uses dada2gpu with GPU support).
// Error models are learned per plate because PCR history affects error rates.

def errExt() { return params.lang == 'python' ? 'pkl' : 'rds' }
def seqExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process DADA2_FILTER_TRIM {
    tag "${meta.id}"
    label 'process_low'
    conda params.lang == 'python' ? "${projectDir}/conda-envs/microscape-python" : "${projectDir}/conda-envs/microscape-r"
    publishDir "${params.outdir}/filtered", mode: 'copy', pattern: "*_filt_stats.tsv", enabled: !params.store_dir

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.filt.fastq.gz"), path("${meta.id}_R2.filt.fastq.gz"), emit: reads
    path("${meta.id}_filt_stats.tsv"), emit: stats

    script:
    if (params.lang == 'python')
    """
    PYTHONPATH=${params.dada2gpu_path}:\$PYTHONPATH \
    dada2_filter_trim.py \
        "${meta.id}" "${r1}" "${r2}" \
        ${params.maxEE} ${params.truncQ} ${params.maxN} \
        ${params.truncLen_fwd} ${params.truncLen_rev} \
        ${task.cpus}
    """
    else
    """
    dada2_filter_trim.R \
        "${meta.id}" "${r1}" "${r2}" \
        ${params.maxEE} ${params.truncQ} ${params.maxN} \
        ${params.truncLen_fwd} ${params.truncLen_rev} \
        ${task.cpus}
    """
}

// Learn error rates from filtered FASTQ files grouped by plate.
process DADA2_LEARN_ERRORS {
    tag "${meta.id}"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/conda-envs/microscape-python" : "${projectDir}/conda-envs/microscape-r"
    publishDir "${params.outdir}/error_models", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/error_models" : null

    input:
    tuple val(meta), path(r1_files), path(r2_files)

    output:
    tuple val(meta.plate), val(meta), path("${meta.id}_errF.${errExt()}"), path("${meta.id}_errR.${errExt()}"), emit: error_models_rds
    path("${meta.id}_error_rates.pdf"), emit: error_plots

    script:
    if (params.lang == 'python')
    """
    PYTHONPATH=${params.dada2gpu_path}:\$PYTHONPATH \
    dada2_learn_errors.py "${meta.id}" ${task.cpus}
    """
    else
    """
    dada2_learn_errors.R "${meta.id}" ${task.cpus}
    """
}

// Per-plate denoising: dereplicate, denoise, merge pairs, build sequence table.
process DADA2_DENOISE {
    tag "${meta.id}"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/conda-envs/microscape-python" : "${projectDir}/conda-envs/microscape-r"
    publishDir "${params.outdir}/seqtabs", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/seqtabs" : null

    input:
    tuple val(meta), path(r1_files), path(r2_files), path(errF), path(errR)

    output:
    tuple val(meta), path("${meta.id}.seqtab.${seqExt()}"), emit: seqtab
    path("${meta.id}.seqtab.tsv"), emit: seqtab_tsv

    script:
    if (params.lang == 'python')
    """
    PYTHONPATH=${params.dada2gpu_path}:\$PYTHONPATH \
    dada2_denoise.py \
        "${meta.id}" "${errF}" "${errR}" \
        ${params.min_overlap} ${task.cpus}
    """
    else
    """
    dada2_denoise.R \
        "${meta.id}" "${errF}" "${errR}" \
        ${params.min_overlap} ${task.cpus}
    """
}

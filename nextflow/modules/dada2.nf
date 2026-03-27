// DADA2 amplicon denoising processes.
//
// Three-step DADA2 workflow:
//   1. DADA2_FILTER_TRIM   — per-sample quality filtering
//   2. DADA2_LEARN_ERRORS  — per-plate error model learning
//   3. DADA2_DENOISE       — per-plate denoising + pair merging → sequence table
//
// Error models are learned per plate because PCR history affects error rates.
// Samples are grouped by plate (first field of sample ID before underscore).
//
// All R logic lives in bin/ scripts; these processes just wire up I/O.

process DADA2_FILTER_TRIM {
    tag "${meta.id}"
    label 'process_low'
    conda "${projectDir}/conda-envs/microscape-dada2"
    publishDir "${params.outdir}/filtered", mode: 'copy', pattern: "*_filt_stats.tsv", enabled: !params.store_dir

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.filt.fastq.gz"), path("${meta.id}_R2.filt.fastq.gz"), emit: reads
    path("${meta.id}_filt_stats.tsv"), emit: stats

    script:
    """
    dada2_filter_trim.R \
        "${meta.id}" "${r1}" "${r2}" \
        ${params.maxEE} ${params.truncQ} ${params.maxN} \
        ${params.truncLen_fwd} ${params.truncLen_rev} \
        ${task.cpus}
    """
}

// Learn error rates from a set of filtered FASTQ files grouped by plate.
// Error models are saved as RDS files for use by DADA2_DENOISE.
process DADA2_LEARN_ERRORS {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/microscape-dada2"
    publishDir "${params.outdir}/error_models", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/error_models" : null

    input:
    tuple val(meta), path(r1_files), path(r2_files)

    output:
    tuple val(meta.plate), val(meta), path("${meta.id}_errF.rds"), path("${meta.id}_errR.rds"), emit: error_models_rds
    path("${meta.id}_error_rates.pdf"), emit: error_plots

    script:
    """
    dada2_learn_errors.R "${meta.id}" ${task.cpus}
    """
}

// Per-plate DADA2 denoising: dereplicate, denoise, merge pairs, build sequence table.
// Uses pre-learned error models from DADA2_LEARN_ERRORS.
process DADA2_DENOISE {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/microscape-dada2"
    publishDir "${params.outdir}/seqtabs", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/seqtabs" : null

    input:
    tuple val(meta), path(r1_files), path(r2_files), path(errF_rds), path(errR_rds)

    output:
    tuple val(meta), path("${meta.id}.seqtab.rds"), emit: seqtab
    path("${meta.id}.seqtab.tsv"), emit: seqtab_tsv

    script:
    """
    dada2_denoise.R \
        "${meta.id}" "${errF_rds}" "${errR_rds}" \
        ${params.min_overlap} ${task.cpus}
    """
}

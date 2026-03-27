// t-SNE ordination of samples and ASVs via Bray-Curtis distances.
//
// Computes Bray-Curtis distance matrices for both samples and ASVs from
// the proportional abundance table, then reduces dimensionality with
// PCA followed by t-SNE for visualization.

process CLUSTER_TSNE {
    tag "tsne"
    label 'process_high'
    conda "${projectDir}/conda-envs/microscape-cluster"
    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    path(seqtab_rds)

    output:
    path("sample_bray_tsne.rds"), emit: sample_tsne
    path("seq_bray_tsne.rds"), emit: seq_tsne
    path("sample_bray_dist.rds"), emit: sample_dist
    path("seq_bray_dist.rds"), emit: seq_dist

    script:
    """
    cluster_tsne.R "${seqtab_rds}" ${task.cpus}
    """
}

// Phylogenetic tree construction.
//
// Aligns unique ASV sequences with DECIPHER and builds a neighbor-joining
// tree. This is optional — only needed for phylogenetic diversity metrics
// (UniFrac, Faith's PD) and phylogeny-aware ordinations.
//
// Can be slow for large ASV sets (>10K sequences). Consider running with
// --store_dir to cache the tree across pipeline runs.

process BUILD_PHYLOGENY {
    tag "phylogeny"
    label 'process_high'
    conda "${projectDir}/conda-envs/microscape-phylo"
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/phylogeny" : null

    input:
    path(seqtab_rds)

    output:
    path("phylo_tree.rds"), emit: tree
    path("phylo_distances.rds"), emit: distances
    path("phylo_seq_map.rds"), emit: seq_map
    path("phylo_alignment.fasta"), emit: alignment

    script:
    """
    build_phylogeny.R "${seqtab_rds}" ${task.cpus}
    """
}

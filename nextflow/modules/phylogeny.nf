// Phylogenetic tree construction.
//
// Aligns unique ASV sequences and builds a neighbor-joining tree.
// R uses DECIPHER; Python uses MAFFT + BioPython/scipy.
// Optional — only needed for UniFrac and phylogeny-aware ordinations.

def phyloExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process BUILD_PHYLOGENY {
    tag "phylogeny"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/conda-envs/microscape-python" : "${projectDir}/conda-envs/microscape-r"
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/phylogeny" : null

    input:
    path(seqtab)

    output:
    path("phylo_tree.${params.lang == 'python' ? 'nwk' : 'rds'}"), emit: tree
    path("phylo_distances.${phyloExt()}"), emit: distances
    path("phylo_seq_map.${phyloExt()}"), emit: seq_map
    path("phylo_alignment.fasta"), emit: alignment

    script:
    if (params.lang == 'python')
    """
    build_phylogeny.py "${seqtab}" ${task.cpus}
    """
    else
    """
    build_phylogeny.R "${seqtab}" ${task.cpus}
    """
}

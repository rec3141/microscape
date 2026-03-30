"""
microscape.phylogeny -- Build MSA and neighbor-joining tree from ASV sequences.

Aligns unique ASV sequences with MAFFT, trims low-coverage positions, computes
a Hamming/gap distance matrix, and builds a neighbor-joining tree.
"""

import os
import tempfile
import shutil
import subprocess
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_fasta(path):
    """Parse a FASTA file into a dict of {id: sequence}."""
    records = {}
    current_id = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if current_id is not None:
            records[current_id] = "".join(parts)
    return records


def _run_mafft(input_fasta, output_fasta, cpus=1):
    """Run MAFFT alignment. Returns True on success, False otherwise."""
    if shutil.which("mafft") is None:
        return False

    cmd = f"mafft --auto --thread {cpus} {input_fasta}"
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=7200
        )
        if result.returncode != 0:
            return False
        with open(output_fasta, "w") as fh:
            fh.write(result.stdout)
        return True
    except (subprocess.TimeoutExpired, Exception):
        return False


def _compute_distance_matrix(trimmed_seqs, n_seqs, trim_len):
    """Compute a pairwise Hamming distance matrix with gap handling."""
    base_to_int = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4, "N": 5}
    int_matrix = np.zeros((n_seqs, trim_len), dtype=np.int8)
    for i, seq in enumerate(trimmed_seqs):
        for j, ch in enumerate(seq):
            int_matrix[i, j] = base_to_int.get(ch, 5)

    dist_matrix = np.zeros((n_seqs, n_seqs), dtype=np.float64)
    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            si = int_matrix[i]
            sj = int_matrix[j]
            both_data = (si < 4) & (sj < 4)
            n_compared = both_data.sum()
            if n_compared == 0:
                dist = 1.0
            else:
                mismatches = ((si != sj) & both_data).sum()
                dist = mismatches / n_compared
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    return dist_matrix


def _build_newick(dist_matrix, id_list, n_seqs):
    """Build a Newick tree string from a distance matrix.

    Tries BioPython's neighbor-joining first, then falls back to scipy
    average-linkage as an NJ approximation.
    """
    newick_str = None

    try:
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
        from Bio import Phylo
        import io

        names = id_list[:]
        matrix_data = []
        for i in range(n_seqs):
            row = [dist_matrix[i, j] for j in range(i + 1)]
            matrix_data.append(row)

        dm = DistanceMatrix(names, matrix_data)
        constructor = DistanceTreeConstructor()
        tree_obj = constructor.nj(dm)

        buf = io.StringIO()
        Phylo.write(tree_obj, buf, "newick")
        newick_str = buf.getvalue().strip()
    except ImportError:
        pass

    if newick_str is None:
        from scipy.cluster.hierarchy import linkage, to_tree
        from scipy.spatial.distance import squareform

        condensed = squareform(dist_matrix)
        Z = linkage(condensed, method="average")

        def _to_newick(node, labels):
            if node.is_leaf():
                return labels[node.id]
            left = _to_newick(node.get_left(), labels)
            right = _to_newick(node.get_right(), labels)
            dl = node.dist / 2
            dr = node.dist / 2
            return f"({left}:{dl:.6f},{right}:{dr:.6f})"

        root = to_tree(Z)
        newick_str = _to_newick(root, id_list) + ";"

    return newick_str


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_phylogeny(sequences, output_newick=None, cpus=1, min_coverage=0.60):
    """Build a multiple sequence alignment and neighbor-joining tree.

    Parameters
    ----------
    sequences : list of str
        DNA sequences (ASVs) to align and build into a tree.
    output_newick : str or None
        If provided, the Newick tree string is also written to this file path.
    cpus : int
        Number of CPU threads to pass to MAFFT.
    min_coverage : float
        Minimum fraction of sequences that must have a non-gap character at
        a given alignment column for it to be retained during trimming.

    Returns
    -------
    dict
        ``"tree_newick"``      -- str, the Newick-format tree.
        ``"distance_matrix"``  -- np.ndarray, pairwise distance matrix.
        ``"alignment"``        -- str, trimmed alignment in FASTA format.
        ``"seq_map"``          -- dict mapping ASV identifiers to original
                                  sequences.
        ``"asv_ids"``          -- list of ASV identifier strings.
    """
    if not sequences:
        raise ValueError("No sequences provided")

    unique_seqs = sorted(set(sequences))
    n_seqs = len(unique_seqs)

    # ------------------------------------------------------------------
    # Assign stable ASV IDs
    # ------------------------------------------------------------------
    seq_map = {}
    for i, seq in enumerate(unique_seqs):
        asv_id = f"ASV_{i + 1:05d}"
        seq_map[asv_id] = seq
    id_list = list(seq_map.keys())

    # ------------------------------------------------------------------
    # Write temp FASTA and run MAFFT
    # ------------------------------------------------------------------
    tmp_dir = tempfile.mkdtemp(prefix="phylo_")
    tmp_input = os.path.join(tmp_dir, "input.fa")
    tmp_aligned = os.path.join(tmp_dir, "aligned.fa")

    try:
        with open(tmp_input, "w") as fh:
            for asv_id in id_list:
                fh.write(f">{asv_id}\n{seq_map[asv_id]}\n")

        mafft_ok = _run_mafft(tmp_input, tmp_aligned, cpus=cpus)

        if not mafft_ok:
            # Fallback: pad sequences to equal length (no true alignment)
            max_len = max(len(s) for s in unique_seqs)
            with open(tmp_aligned, "w") as fh:
                for asv_id in id_list:
                    padded = seq_map[asv_id].ljust(max_len, "-")
                    fh.write(f">{asv_id}\n{padded}\n")

        # --------------------------------------------------------------
        # Parse aligned FASTA
        # --------------------------------------------------------------
        aligned = _parse_fasta(tmp_aligned)
        aligned_seqs = [aligned[aid].upper() for aid in id_list]
        aln_len = max(len(s) for s in aligned_seqs) if aligned_seqs else 0
        aligned_seqs = [s.ljust(aln_len, "-") for s in aligned_seqs]

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Trim alignment (positions where < min_coverage of seqs have data)
    # ------------------------------------------------------------------
    aln_matrix = np.array([list(s) for s in aligned_seqs])
    coverage = np.mean(aln_matrix != "-", axis=0)
    keep_cols = coverage >= min_coverage
    n_kept = keep_cols.sum()

    if n_kept == 0:
        raise RuntimeError(
            "All alignment positions were trimmed -- no informative data remaining"
        )

    aln_matrix = aln_matrix[:, keep_cols]
    trimmed_seqs = ["".join(row) for row in aln_matrix]
    trim_len = len(trimmed_seqs[0])

    # ------------------------------------------------------------------
    # Compute distance matrix
    # ------------------------------------------------------------------
    dist_matrix = _compute_distance_matrix(trimmed_seqs, n_seqs, trim_len)

    # ------------------------------------------------------------------
    # Build tree
    # ------------------------------------------------------------------
    newick_str = _build_newick(dist_matrix, id_list, n_seqs)

    # ------------------------------------------------------------------
    # Build alignment FASTA string
    # ------------------------------------------------------------------
    alignment_parts = []
    for aid, seq in zip(id_list, trimmed_seqs):
        alignment_parts.append(f">{aid}\n{seq}")
    alignment_str = "\n".join(alignment_parts) + "\n"

    # ------------------------------------------------------------------
    # Optionally write Newick to file
    # ------------------------------------------------------------------
    if output_newick is not None:
        with open(output_newick, "w") as fh:
            fh.write(newick_str + "\n")

    return {
        "tree_newick": newick_str,
        "distance_matrix": dist_matrix,
        "alignment": alignment_str,
        "seq_map": seq_map,
        "asv_ids": id_list,
    }

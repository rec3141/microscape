#!/usr/bin/env python3
#
# build_phylogeny.py — Build MSA and neighbor-joining tree from ASV sequences
#
# Aligns unique ASV sequences with MAFFT, trims low-coverage positions,
# computes a Hamming/gap distance matrix, and builds a neighbor-joining tree.
#
# Usage:
#   build_phylogeny.py <seqtab.pkl> <cpus>
#
# Outputs:
#   phylo_tree.nwk          Newick tree file
#   phylo_tree.pkl          tree object (Bio.Phylo)
#   phylo_distances.pkl     distance matrix (DataFrame)
#   phylo_seq_map.pkl       ASV_id -> sequence mapping
#   phylo_alignment.fasta   multiple sequence alignment

import sys
import os
import pickle
import time
import tempfile
import subprocess
import shutil
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if len(sys.argv) < 3:
    print("Usage: build_phylogeny.py <seqtab.pkl> <cpus>", file=sys.stderr)
    sys.exit(1)

input_path = sys.argv[1]
cpus = int(sys.argv[2])

t0 = time.time()

# ---------------------------------------------------------------------------
# Load sequences
# ---------------------------------------------------------------------------
print(f"[INFO] Loading sequences from {input_path}")
with open(input_path, "rb") as fh:
    data = pickle.load(fh)

if isinstance(data, pd.DataFrame):
    if "sequence" in data.columns:
        unique_seqs = sorted(data["sequence"].unique().tolist())
    else:
        print("[ERROR] DataFrame has no 'sequence' column", file=sys.stderr)
        sys.exit(1)
elif isinstance(data, list):
    unique_seqs = sorted(set(data))
else:
    print("[ERROR] Unsupported input type: " + str(type(data)), file=sys.stderr)
    sys.exit(1)

n_seqs = len(unique_seqs)
print(f"[INFO] {n_seqs} unique sequences to align")

if n_seqs == 0:
    print("[ERROR] No sequences found in input", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Assign stable ASV IDs
# ---------------------------------------------------------------------------
seq_map = {}
for i, seq in enumerate(unique_seqs):
    asv_id = f"ASV_{i+1:05d}"
    seq_map[asv_id] = seq

id_list = list(seq_map.keys())
print(f"[INFO] Assigned {len(id_list)} ASV identifiers")

# ---------------------------------------------------------------------------
# Write temp FASTA for MAFFT
# ---------------------------------------------------------------------------
tmp_dir = tempfile.mkdtemp(prefix="phylo_")
tmp_input = os.path.join(tmp_dir, "input.fa")
tmp_aligned = os.path.join(tmp_dir, "aligned.fa")

with open(tmp_input, "w") as fh:
    for asv_id in id_list:
        fh.write(f">{asv_id}\n{seq_map[asv_id]}\n")

# ---------------------------------------------------------------------------
# Run MAFFT
# ---------------------------------------------------------------------------
mafft_available = shutil.which("mafft") is not None

if mafft_available:
    print(f"[INFO] Running MAFFT with --auto --thread {cpus} ...")
    cmd = f"mafft --auto --thread {cpus} {tmp_input}"
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=7200
        )
        if result.returncode != 0:
            print(f"[ERROR] MAFFT failed: {result.stderr[:500]}", file=sys.stderr)
            mafft_available = False
        else:
            with open(tmp_aligned, "w") as fh:
                fh.write(result.stdout)
            print(f"[INFO] MAFFT alignment complete")
    except subprocess.TimeoutExpired:
        print("[ERROR] MAFFT timed out after 7200s", file=sys.stderr)
        mafft_available = False
    except Exception as e:
        print(f"[ERROR] MAFFT error: {e}", file=sys.stderr)
        mafft_available = False

if not mafft_available:
    # Fallback: no alignment, just pad sequences to equal length
    print("[INFO] MAFFT not available — falling back to unaligned padded sequences")
    max_len = max(len(s) for s in unique_seqs)
    with open(tmp_aligned, "w") as fh:
        for asv_id in id_list:
            padded = seq_map[asv_id].ljust(max_len, "-")
            fh.write(f">{asv_id}\n{padded}\n")
    print("[INFO] Padded sequences written (no true alignment)")

# ---------------------------------------------------------------------------
# Parse aligned FASTA
# ---------------------------------------------------------------------------
def parse_fasta(path):
    """Return dict {id: sequence}."""
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

aligned = parse_fasta(tmp_aligned)
print(f"[INFO] Parsed alignment: {len(aligned)} sequences")

# Ensure ordering matches id_list
aligned_seqs = [aligned[aid].upper() for aid in id_list]
aln_len = max(len(s) for s in aligned_seqs) if aligned_seqs else 0
# Pad any short ones (shouldn't happen, but be safe)
aligned_seqs = [s.ljust(aln_len, "-") for s in aligned_seqs]

print(f"[INFO] Alignment length: {aln_len} positions")

# ---------------------------------------------------------------------------
# Trim alignment ends (positions where <60% of sequences have data)
# ---------------------------------------------------------------------------
MIN_COVERAGE = 0.60

aln_matrix = np.array([list(s) for s in aligned_seqs])  # (n_seqs, aln_len)
coverage = np.mean(aln_matrix != "-", axis=0)            # fraction non-gap per column

keep_cols = coverage >= MIN_COVERAGE
n_kept = keep_cols.sum()
print(f"[INFO] Trimming: keeping {n_kept}/{aln_len} positions "
      f"(>={MIN_COVERAGE*100:.0f}% coverage)")

if n_kept == 0:
    print("[ERROR] All alignment positions trimmed — no data remaining", file=sys.stderr)
    sys.exit(1)

aln_matrix = aln_matrix[:, keep_cols]
trimmed_seqs = ["".join(row) for row in aln_matrix]
trim_len = len(trimmed_seqs[0])

# ---------------------------------------------------------------------------
# Compute distance matrix (Hamming + gap handling)
# ---------------------------------------------------------------------------
print("[INFO] Computing pairwise distance matrix ...")

dist_matrix = np.zeros((n_seqs, n_seqs), dtype=np.float64)

# Convert to numeric for fast comparison
base_to_int = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4, "N": 5}
int_matrix = np.zeros((n_seqs, trim_len), dtype=np.int8)
for i, seq in enumerate(trimmed_seqs):
    for j, ch in enumerate(seq):
        int_matrix[i, j] = base_to_int.get(ch, 5)

for i in range(n_seqs):
    for j in range(i + 1, n_seqs):
        si = int_matrix[i]
        sj = int_matrix[j]
        # Only compare positions where both have data (not gap)
        both_data = (si < 4) & (sj < 4)
        n_compared = both_data.sum()
        if n_compared == 0:
            dist = 1.0
        else:
            mismatches = ((si != sj) & both_data).sum()
            dist = mismatches / n_compared
        dist_matrix[i, j] = dist
        dist_matrix[j, i] = dist

dist_df = pd.DataFrame(dist_matrix, index=id_list, columns=id_list)
print(f"[INFO] Distance matrix computed ({n_seqs}x{n_seqs})")

# ---------------------------------------------------------------------------
# Build neighbor-joining tree
# ---------------------------------------------------------------------------
print("[INFO] Building neighbor-joining tree ...")

tree_obj = None
newick_str = None

try:
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from Bio import Phylo
    import io

    # Bio.Phylo wants a lower-triangular list of lists
    names = id_list[:]
    matrix_data = []
    for i in range(n_seqs):
        row = [dist_matrix[i, j] for j in range(i + 1)]
        matrix_data.append(row)

    dm = DistanceMatrix(names, matrix_data)
    constructor = DistanceTreeConstructor()
    tree_obj = constructor.nj(dm)

    # Convert to Newick string
    buf = io.StringIO()
    Phylo.write(tree_obj, buf, "newick")
    newick_str = buf.getvalue().strip()
    print("[INFO] NJ tree built with BioPython")

except ImportError:
    print("[INFO] BioPython not available — building NJ tree with scipy")
    try:
        from scipy.cluster.hierarchy import linkage, to_tree
        from scipy.spatial.distance import squareform

        condensed = squareform(dist_matrix)
        # Use average linkage as NJ approximation
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
        print("[INFO] NJ tree built with scipy (average linkage approximation)")
    except ImportError:
        print("[ERROR] Neither BioPython nor scipy available for tree building",
              file=sys.stderr)
        sys.exit(1)

# Try to parse newick into Bio.Phylo tree object for pickle
if tree_obj is None and newick_str is not None:
    try:
        from Bio import Phylo
        import io
        tree_obj = Phylo.read(io.StringIO(newick_str), "newick")
    except Exception:
        tree_obj = newick_str  # Store as string if Phylo unavailable

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------

# 1. Newick tree
with open("phylo_tree.nwk", "w") as fh:
    fh.write(newick_str + "\n")
print("[INFO] Saved phylo_tree.nwk")

# 2. Tree object pickle
with open("phylo_tree.pkl", "wb") as fh:
    pickle.dump(tree_obj, fh, protocol=4)
print("[INFO] Saved phylo_tree.pkl")

# 3. Distance matrix
with open("phylo_distances.pkl", "wb") as fh:
    pickle.dump(dist_df, fh, protocol=4)
print("[INFO] Saved phylo_distances.pkl")

# 4. Sequence map
with open("phylo_seq_map.pkl", "wb") as fh:
    pickle.dump(seq_map, fh, protocol=4)
print("[INFO] Saved phylo_seq_map.pkl")

# 5. Alignment FASTA (trimmed)
with open("phylo_alignment.fasta", "w") as fh:
    for aid, seq in zip(id_list, trimmed_seqs):
        fh.write(f">{aid}\n{seq}\n")
print("[INFO] Saved phylo_alignment.fasta")

# Cleanup temp files
shutil.rmtree(tmp_dir, ignore_errors=True)

elapsed = time.time() - t0
print(f"[INFO] Phylogeny pipeline complete in {elapsed:.1f}s")

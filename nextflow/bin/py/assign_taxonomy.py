#!/usr/bin/env python3
#
# assign_taxonomy.py — Taxonomy assignment using dada2's C-level classifier
#
# Calls the same C code as R's dada2::assignTaxonomy() via ctypes, ensuring
# identical results. The C library uses 8-mer naive Bayesian classification
# with 100 bootstrap iterations (Wang et al. 2007).
#
# The heavy lifting (kmer table construction, log-probability scoring,
# bootstrap resampling) runs in compiled C with OpenMP parallelism.
# Python just handles I/O and the reference database parsing.
#
# Usage:
#   assign_taxonomy.py <seqtab.pkl> <ref_db.fasta> <db_name> <cpus> [tax_levels]
#
# Outputs:
#   <db_name>_taxonomy.pkl     dict with 'tax' and 'boot' DataFrames
#   <db_name>_bootstrap.pkl    bootstrap DataFrame
#   <db_name>_taxonomy.tsv     tab-delimited taxonomy

import sys
import os
import pickle
import gzip
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if len(sys.argv) < 5:
    print("Usage: assign_taxonomy.py <seqtab.pkl> <ref_db> <db_name> <cpus> [tax_levels]",
          file=sys.stderr)
    sys.exit(1)

input_path = sys.argv[1]
ref_path   = sys.argv[2]
db_name    = sys.argv[3]
cpus       = int(sys.argv[4])
tax_levels = sys.argv[5].split(",") if len(sys.argv) > 5 and sys.argv[5] != "null" else None

os.environ["OMP_NUM_THREADS"] = str(cpus)

# ---------------------------------------------------------------------------
# Load query sequences
# ---------------------------------------------------------------------------
print(f"[INFO] Loading query data from {input_path}")
data = pickle.load(open(input_path, "rb"))
if isinstance(data, pd.DataFrame):
    seqs = sorted(data["sequence"].unique())
elif isinstance(data, list):
    seqs = data
else:
    seqs = list(data)
print(f"[INFO] {len(seqs)} unique query sequences")

# ---------------------------------------------------------------------------
# Parse reference database
# ---------------------------------------------------------------------------
print(f"[INFO] Parsing reference database: {ref_path}")
ref_seqs = []
ref_tax_tokens = []
current_seq_parts = []
current_tax = None

_open = gzip.open if ref_path.endswith(".gz") else open
with _open(ref_path, "rt") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if current_tax is not None:
                ref_seqs.append("".join(current_seq_parts).upper().replace("U", "T"))
                ref_tax_tokens.append(current_tax)
            tokens = [t.strip() for t in line[1:].split(";")]
            tokens = [t for t in tokens if t]
            current_tax = tokens
            current_seq_parts = []
        else:
            current_seq_parts.append(line.strip())
    if current_tax is not None:
        ref_seqs.append("".join(current_seq_parts).upper().replace("U", "T"))
        ref_tax_tokens.append(current_tax)

nref = len(ref_seqs)
print(f"[INFO] {nref} reference sequences")

# ---------------------------------------------------------------------------
# Determine taxonomy levels and build genus matrix
#
# The C code needs:
#   - ref_to_genus: mapping each reference to a unique genus ID (0-indexed)
#   - genusmat: genus-to-taxonomy-level assignment matrix (ngenus x nlevel)
#     where each cell is an integer ID for that taxon at that level
# ---------------------------------------------------------------------------
max_depth = max(len(t) for t in ref_tax_tokens)
if tax_levels is None:
    nlevel = max_depth
    tax_levels = [f"Level_{i+1}" for i in range(nlevel)]
else:
    nlevel = len(tax_levels)

# Pad all taxonomy tokens to uniform depth
for i in range(len(ref_tax_tokens)):
    while len(ref_tax_tokens[i]) < nlevel:
        ref_tax_tokens[i].append("")

# Build the genus map: unique taxonomy tuples → genus ID
# "Genus" here means the finest taxonomic level (last column of genusmat)
genus_map = {}  # tuple of taxonomy tokens → genus_id
genus_list = []  # list of taxonomy token tuples

ref_to_genus = np.zeros(nref, dtype=np.int32)
for i in range(nref):
    key = tuple(ref_tax_tokens[i][:nlevel])
    if key not in genus_map:
        genus_map[key] = len(genus_list)
        genus_list.append(key)
    ref_to_genus[i] = genus_map[key]

ngenus = len(genus_list)
print(f"[INFO] {ngenus} unique genus-level groups across {nlevel} levels")

# Build genusmat: encode each taxon string as an integer
# At each level, map unique strings to integer IDs
level_maps = [{} for _ in range(nlevel)]
genusmat = np.zeros((ngenus, nlevel), dtype=np.int32)

for g, tax_tuple in enumerate(genus_list):
    for lev in range(nlevel):
        taxon = tax_tuple[lev] if lev < len(tax_tuple) else ""
        if taxon not in level_maps[lev]:
            level_maps[lev][taxon] = len(level_maps[lev])
        genusmat[g, lev] = level_maps[lev][taxon]

# Build reverse maps for decoding results
level_reverse = []
for lev in range(nlevel):
    rev = {v: k for k, v in level_maps[lev].items()}
    level_reverse.append(rev)

# ---------------------------------------------------------------------------
# Call C library for taxonomy assignment
# ---------------------------------------------------------------------------
print(f"[INFO] Running C-level taxonomy classifier (k=8, {cpus} threads)...")

# Import after setting OMP_NUM_THREADS
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "..", "dada2_gpu"))
try:
    from dada2gpu._cdada import run_taxonomy
except ImportError:
    # Try with PYTHONPATH
    from dada2gpu._cdada import run_taxonomy

result = run_taxonomy(
    seqs=seqs,
    refs=ref_seqs,
    ref_to_genus=ref_to_genus,
    genusmat=genusmat,
    ngenus=ngenus,
    nlevel=nlevel,
    verbose=True
)

# ---------------------------------------------------------------------------
# Decode results into taxonomy strings
# ---------------------------------------------------------------------------
rval = result["rval"]    # (nseq,) 1-indexed genus ID, 0=NA
rboot = result["rboot"]  # (nseq, nlevel) bootstrap counts out of 100

tax_dict = {}
boot_dict = {}

for i, seq in enumerate(seqs):
    g = rval[i] - 1  # to 0-indexed
    if g < 0 or g >= ngenus:
        # No assignment
        tax_dict[seq] = {tax_levels[lev]: "" for lev in range(nlevel)}
        boot_dict[seq] = {tax_levels[lev]: 0 for lev in range(nlevel)}
    else:
        genus_tax = genus_list[g]
        tax_dict[seq] = {tax_levels[lev]: genus_tax[lev] if lev < len(genus_tax) else ""
                         for lev in range(nlevel)}
        boot_dict[seq] = {tax_levels[lev]: int(rboot[i, lev])
                          for lev in range(nlevel)}

tax_df = pd.DataFrame.from_dict(tax_dict, orient="index")
boot_df = pd.DataFrame.from_dict(boot_dict, orient="index")

# Summary
for col in tax_df.columns:
    n_classified = (tax_df[col] != "").sum()
    pct = round(n_classified / len(tax_df) * 100, 1) if len(tax_df) > 0 else 0
    print(f"[INFO] {col}: {n_classified}/{len(tax_df)} ({pct}%) classified")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
with open(f"{db_name}_taxonomy.pkl", "wb") as f:
    pickle.dump({"tax": tax_df, "boot": boot_df}, f)

with open(f"{db_name}_bootstrap.pkl", "wb") as f:
    pickle.dump(boot_df, f)

tax_out = tax_df.copy()
tax_out.insert(0, "sequence", tax_out.index)
tax_out.to_csv(f"{db_name}_taxonomy.tsv", sep="\t", index=False)

print(f"[INFO] Taxonomy assignment complete for {db_name}")

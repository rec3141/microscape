#!/bin/bash -ue
PYTHONPATH=/data/microscape_v2/nextflow/../dada2_gpu:${PYTHONPATH:-}     assign_taxonomy.py         "seqtab_final.pkl" "ref_dada2_silva_v138-2_train_set_uniq.fasta.gz" "silva"         8 "Domain,Phylum,Class,Order,Family,Genus"

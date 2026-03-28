#!/bin/bash -ue
PYTHONPATH=/data/microscape_v2/nextflow/../dada2_gpu:${PYTHONPATH:-}     remove_chimeras.py "seqtab_merged.pkl" 8

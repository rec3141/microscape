#!/bin/bash -ue
PYTHONPATH=/data/microscape_v2/nextflow/../dada2_gpu:${PYTHONPATH:-}     dada2_learn_errors.py "16S-Plate1Pool-Rep1" 8

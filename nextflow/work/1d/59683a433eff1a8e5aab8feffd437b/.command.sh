#!/bin/bash -ue
PYTHONPATH=/data/microscape_v2/nextflow/../dada2_gpu:${PYTHONPATH:-}     dada2_denoise.py         "16S-Plate1Pool-Rep1" "16S-Plate1Pool-Rep1_errF.pkl" "16S-Plate1Pool-Rep1_errR.pkl"         10 8

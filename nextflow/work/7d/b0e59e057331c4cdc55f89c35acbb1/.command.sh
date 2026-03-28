#!/bin/bash -ue
PYTHONPATH=/data/microscape_v2/nextflow/../dada2_gpu:${PYTHONPATH:-}     dada2_filter_trim.py         "16S-Plate1Pool-Rep1_G11" "16S-Plate1Pool-Rep1_G11_R1.fastq.gz" "16S-Plate1Pool-Rep1_G11_R2.fastq.gz"         2 11 0         0 0         2

#!/bin/bash

source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
conda activate skipper
snakemake -kps Skipper.py -j 30 -w 30

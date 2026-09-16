#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=Test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
conda activate skipper
ulimit -n 65536

snakemake -ns Skipper.py -j 1


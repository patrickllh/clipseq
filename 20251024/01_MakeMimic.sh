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

# get random 60% reads to make mimic replicate
seqkit sample -p 0.6 input/N288_L1_1.fq.gz -o input/mimic_N288_L1_1.fq.gz
seqkit sample -p 0.6 input/N289_L1_1.fq.gz -o input/mimic_N289_L1_1.fq.gz
seqkit sample -p 0.6 input/N290_L1_1.fq.gz -o input/mimic_N290_L1_1.fq.gz
seqkit sample -p 0.6 input/N291_L1_1.fq.gz -o input/mimic_N291_L1_1.fq.gz


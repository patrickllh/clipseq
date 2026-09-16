#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=02c_STAR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate skipper
conda list

mkdir -p output

STAR --runMode genomeGenerate --runThreadN 52 \
--genomeDir output/genome_ref \
--genomeFastaFiles input/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--sjdbGTFtagExonParentTranscript input/gencode.v38.annotation.gff3.gz

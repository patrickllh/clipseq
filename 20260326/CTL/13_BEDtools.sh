#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=13_BEDtools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate merge_peaks
conda list

mkdir -p output

bedtools intersect \
-a output/SNRNP70-Con-rep1_L1_R1.peaks.normed.compressed.bed \
-b output/SNRNP70-Con-rep2_L2_R1.peaks.normed.compressed.bed \
| sort -k1,1 -k2,2n \
| bedtools merge -d 0 \
| bedtools sort -i - -g output/genome_ref/chrNameLength.txt \
> output/SNRNP70-Con.peaks.intersect.bed

cat output/mVenus-Con-rep5_L1_R1.peaks.bed output/mVenus-Con-rep6_L2_R1.peaks.bed \
| sort -k1,1 -k2,2n \
| bedtools merge \
> output/mVenus-Con.peaks.blacklist.bed

bedtools intersect -v \
-a output/SNRNP70-Con.peaks.intersect.bed \
-b output/mVenus-Con.peaks.blacklist.bed \
> output/SNRNP70-Con.peaks.consensus.bed

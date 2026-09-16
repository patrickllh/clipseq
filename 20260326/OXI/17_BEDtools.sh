#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=17_BEDtools
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

bedgraphs=(
SNRNP70-SA-rep1_L1_R1.cpm.bw.bedGraph
SNRNP70-SA-rep2_L2_R1.cpm.bw.bedGraph
mVenus-SA-rep5_L1_R1.cpm.bw.bedGraph
mVenus-SA-rep6_L2_R1.cpm.bw.bedGraph
)

for bg in "${bedgraphs[@]}"; do
    FILE="output/${bg}"
    bedtools intersect -a "${FILE}" -b output/SNRNP70-SA.peaks.consensus.bed > "output/${bg}.inside.bedGraph"
    bedtools complement -i output/SNRNP70-SA.peaks.consensus.bed -g output/genome_ref/chrNameLength.txt > "output/${bg}.outside.bed"
    awk '{print $1, $2, $3, 0}' OFS='\t' "output/${bg}.outside.bed" > "output/${bg}.zeros.bedGraph"
    cat "output/${bg}.inside.bedGraph" "output/${bg}.zeros.bedGraph" | sort -k1,1 -k2,2n > "output/${bg}.final.bedGraph"
done

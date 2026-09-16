#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=16_UCSC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate ucsc_tools
conda list

mkdir -p output

bigwigs=(
SNRNP70-SA-rep1_L1_R1.cpm.bw
SNRNP70-SA-rep2_L2_R1.cpm.bw
mVenus-SA-rep5_L1_R1.cpm.bw
mVenus-SA-rep6_L2_R1.cpm.bw
)

for bw in "${bigwigs[@]}"; do
    FILE="output/${bw}"
    bigWigToBedGraph "${FILE}" "output/${bw}.bedGraph"
done

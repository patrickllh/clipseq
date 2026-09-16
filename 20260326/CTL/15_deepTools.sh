#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=15_deepTools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate deeptools
conda list

mkdir -p output

multiBigwigSummary BED-file --bwfiles \
output/SNRNP70-Con-rep1_L1_R1.cpm.bw \
output/SNRNP70-Con-rep2_L2_R1.cpm.bw \
output/mVenus-Con-rep5_L1_R1.cpm.bw \
output/mVenus-Con-rep6_L2_R1.cpm.bw \
--BED output/SNRNP70-Con.peaks.consensus.bed \
--outFileName output/SNRNP70-Con.peaks.consensus.cpm.average.npz \
--outRawCounts output/SNRNP70-Con.peaks.consensus.cpm.average.tab

# i set pseudocount: +1 to both IP and IN
awk 'BEGIN{
  OFS="\t";
  print "chr","start","end","fc1","fc2","mean_fc","log2_mean_fc"
}
NR>1 {
  ip1 = $4 + 1
  ip2 = $5 + 1
  in1 = $6 + 1
  in2 = $7 + 1

  fc1 = ip1 / in1
  fc2 = ip2 / in2
  mean_fc = (fc1 + fc2) / 2
  log2fc = log(mean_fc)/log(2)

  print $1,$2,$3,fc1,fc2,mean_fc,log2fc
}' output/SNRNP70-Con.peaks.consensus.cpm.average.tab \
> output/SNRNP70-Con.peaks.consensus.foldchange.tab

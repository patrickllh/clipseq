#!/bin/bash
#SBATCH --partition=cluster_long
#SBATCH --job-name=multiBigwigSummary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
conda activate deeptools

multiBigwigSummary bins \
  -b \
  ~/clipseq/20260421/CTL/output/SNRNP70-Con-rep1_L1_R1.cpm.bw \
  ~/clipseq/20260421/CTL/output/SNRNP70-Con-rep2_L2_R1.cpm.bw \
  ~/clipseq/20260421/CTL/output/mVenus-Con-rep5_L1_R1.cpm.bw \
  ~/clipseq/20260421/CTL/output/mVenus-Con-rep6_L2_R1.cpm.bw \
  ~/clipseq/20260421/OSM/output/SNRNP70-Sorbitol-rep1_L1_R1.cpm.bw \
  ~/clipseq/20260421/OSM/output/SNRNP70-Sorbitol-rep2_L2_R1.cpm.bw \
  ~/clipseq/20260421/OSM/output/mVenus-Sorbitol-rep5_L1_R1.cpm.bw \
  ~/clipseq/20260421/OSM/output/mVenus-Sorbitol-rep6_L2_R1.cpm.bw \
  -o SNRNP70_mVenus_CTL_OSM.bw_summary.npz \
  --labels \
  SNRNP70-Con-rep1 SNRNP70-Con-rep2 mVenus-Con-rep5 mVenus-Con-rep6 \
  SNRNP70-Sorbitol-rep1 SNRNP70-Sorbitol-rep2 mVenus-Sorbitol-rep5 mVenus-Sorbitol-rep6 \
  -p 48

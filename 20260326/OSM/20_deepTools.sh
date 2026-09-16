#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=20_deepTools
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

computeMatrix scale-regions --metagene -p max --verbose \
-a 1000 -b 1000 -R output/mrna.gtf \
-S \
output/SNRNP70-Sorbitol-rep1_L1_R1.cpm.peaks.consensus.bw \
output/SNRNP70-Sorbitol-rep2_L2_R1.cpm.peaks.consensus.bw \
output/mVenus-Sorbitol-rep5_L1_R1.cpm.peaks.consensus.bw \
output/mVenus-Sorbitol-rep6_L2_R1.cpm.peaks.consensus.bw \
-o output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.mat.gz

plotHeatmap --linesAtTickMark --colorMap Blues --heatmapHeight 10 --zMax 0.5 \
--samplesLabel OSM_SmV_1 OSM_SmV_2 OSM_mV_1 OSM_mV_2 \
-m output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.mat.gz \
-o output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.heat.pdf \
--outFileSortedRegions output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.heat.bed \
--outFileNameMatrix output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.heat.mat.gz

computeMatrix scale-regions --metagene -p max --verbose \
-a 1000 -b 1000 --skipZeros -R output/mrna.gtf \
-S \
output/SNRNP70-Sorbitol-rep1_L1_R1.cpm.peaks.consensus.bw \
output/SNRNP70-Sorbitol-rep2_L2_R1.cpm.peaks.consensus.bw \
output/mVenus-Sorbitol-rep5_L1_R1.cpm.peaks.consensus.bw \
output/mVenus-Sorbitol-rep6_L2_R1.cpm.peaks.consensus.bw \
-o output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.skipzeros.mat.gz

plotHeatmap --linesAtTickMark --colorMap Blues --heatmapHeight 10 --zMax 0.5 \
--samplesLabel OSM_SmV_1 OSM_SmV_2 OSM_mV_1 OSM_mV_2 \
-m output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.skipzeros.mat.gz \
-o output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.skipzeros.heat.pdf \
--outFileSortedRegions output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.skipzeros.heat.bed \
--outFileNameMatrix output/mRNA.SNRNP70-Sorbitol.cpm.peaks.consensus.metagene.skipzeros.heat.mat.gz

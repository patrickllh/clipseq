#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=01b_deepTools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate deeptools

computeMatrix scale-regions --metagene -p max --verbose \
--skipZeros -R ~/clipseq/20260421/extract-transcript-regions/gtf_extract_transcript_data/transcript_regions.bed \
-S \
~/clipseq/20260421/OSM/output/SNRNP70-Sorbitol-rep1_L1_R1.crosslinks.cpm.bw \
~/clipseq/20260421/OSM/output/mVenus-Sorbitol-rep5_L1_R1.crosslinks.cpm.bw \
-o SNRNP70-Sorbitol-rep1.crosslinks.mrna.skipzeros.mat.gz

mrna=$(zcat SNRNP70-Sorbitol-rep1.crosslinks.mrna.skipzeros.mat.gz | wc -l)

plotHeatmap --colorMap Blues --heatmapHeight 10 --legendLocation none \
--regionsLabel "${mrna} mRNAs" \
--startLabel "TSS" \
--endLabel "TES" \
--yAxisLabel "mean RPM" \
--xAxisLabel "Scaled regions" \
--samplesLabel OSM_SmV_1 OSM_mV_1 \
--sortUsingSamples 1 \
-m SNRNP70-Sorbitol-rep1.crosslinks.mrna.skipzeros.mat.gz \
-o SNRNP70-Sorbitol-rep1.crosslinks.mrna.skipzeros.heat.pdf \
--outFileSortedRegions SNRNP70-Sorbitol-rep1.crosslinks.mrna.skipzeros.heat.bed

computeMatrix scale-regions --metagene -p max --verbose \
--skipZeros -R ~/clipseq/20260421/extract-transcript-regions/gtf_extract_transcript_data/transcript_regions.bed \
-S \
~/clipseq/20260421/OSM/output/SNRNP70-Sorbitol-rep2_L2_R1.crosslinks.cpm.bw \
~/clipseq/20260421/OSM/output/mVenus-Sorbitol-rep6_L2_R1.crosslinks.cpm.bw \
-o SNRNP70-Sorbitol-rep2.crosslinks.mrna.skipzeros.mat.gz

mrna=$(zcat SNRNP70-Sorbitol-rep2.crosslinks.mrna.skipzeros.mat.gz | wc -l)

plotHeatmap --colorMap Blues --heatmapHeight 10 --legendLocation none \
--regionsLabel "${mrna} mRNAs" \
--startLabel "TSS" \
--endLabel "TES" \
--yAxisLabel "mean RPM" \
--xAxisLabel "Scaled regions" \
--samplesLabel OSM_SmV_2 OSM_mV_2 \
--sortUsingSamples 1 \
-m SNRNP70-Sorbitol-rep2.crosslinks.mrna.skipzeros.mat.gz \
-o SNRNP70-Sorbitol-rep2.crosslinks.mrna.skipzeros.heat.pdf \
--outFileSortedRegions SNRNP70-Sorbitol-rep2.mrna.crosslinks.skipzeros.heat.bed

#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=01_deepTools
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

computeMatrix scale-regions --metagene -p max --verbose \
--skipZeros -R ~/extract-transcript-regions/gtf_extract_transcript_data/transcript_regions.bed \
-S \
~/clipseq/20260421/CTL/output/SNRNP70-Con-rep1_L1_R1.cpm.bw \
~/clipseq/20260421/CTL/output/mVenus-Con-rep5_L1_R1.cpm.bw \
-o SNRNP70-Con-rep1.mrna.skipzeros.mat.gz

mrna=$(zcat SNRNP70-Con-rep1.mrna.skipzeros.mat.gz | wc -l)

plotHeatmap --colorMap Blues --heatmapHeight 10 --legendLocation none \
--regionsLabel "${mrna} mRNAs" \
--startLabel "TSS" \
--endLabel "TES" \
--yAxisLabel "mean RPM" \
--xAxisLabel "Scaled regions" \
--samplesLabel CTL_SmV_1 CTL_mV_1 \
--sortUsingSamples 1 \
-m SNRNP70-Con-rep1.mrna.skipzeros.mat.gz \
-o SNRNP70-Con-rep1.mrna.skipzeros.heat.pdf \
--outFileSortedRegions SNRNP70-Con-rep1.mrna.skipzeros.heat.bed

computeMatrix scale-regions --metagene -p max --verbose \
--skipZeros -R ~/extract-transcript-regions/gtf_extract_transcript_data/transcript_regions.bed \
-S \
~/clipseq/20260421/CTL/output/SNRNP70-Con-rep2_L2_R1.cpm.bw \
~/clipseq/20260421/CTL/output/mVenus-Con-rep6_L2_R1.cpm.bw \
-o SNRNP70-Con-rep2.mrna.skipzeros.mat.gz

mrna=$(zcat SNRNP70-Con-rep2.mrna.skipzeros.mat.gz | wc -l)

plotHeatmap --colorMap Blues --heatmapHeight 10 --legendLocation none \
--regionsLabel "${mrna} mRNAs" \
--startLabel "TSS" \
--endLabel "TES" \
--yAxisLabel "mean RPM" \
--xAxisLabel "Scaled regions" \
--samplesLabel CTL_SmV_2 CTL_mV_2 \
--sortUsingSamples 1 \
-m SNRNP70-Con-rep2.mrna.skipzeros.mat.gz \
-o SNRNP70-Con-rep2.mrna.skipzeros.heat.pdf \
--outFileSortedRegions SNRNP70-Con-rep2.mrna.skipzeros.heat.bed

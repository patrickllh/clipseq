#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=12_merge_peaks
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

# note: please make each replicate correspond to the same order per array
ip_reads=(
SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.sort.dedup.bam
SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.sort.dedup.bam
)

in_reads=(
mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.sort.dedup.bam
mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.sort.dedup.bam
)

ip_peaks=(
SNRNP70-Sorbitol-rep1_L1_R1.peaks.bed
SNRNP70-Sorbitol-rep2_L2_R1.peaks.bed
)

ip_readnum=(
SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.sort.dedup.bam.readnum
SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.sort.dedup.bam.readnum
)

in_readnum=(
mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.sort.dedup.bam.readnum
mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.sort.dedup.bam.readnum
)

for i in "${!ip_reads[@]}"; do
    sample=$(basename "${ip_peaks[$i]}" .bed)
    perl ~/merge_peaks/bin/perl/overlap_peakfi_with_bam.pl \
        output/${ip_reads[$i]} \
        output/${in_reads[$i]} \
        output/${ip_peaks[$i]} \
        output/${ip_readnum[$i]} \
        output/${in_readnum[$i]} \
        "output/${sample}.normed.bed"
    perl ~/merge_peaks/bin/perl/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
        "output/${sample}.normed.bed.full" \
        "output/${sample}.normed.compressed.bed" \
        "output/${sample}.normed.compressed.bed.full"
done

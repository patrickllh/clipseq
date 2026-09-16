#!/bin/bash
#SBATCH --partition=cluster_long
#SBATCH --job-name=bamCoverage.sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

# get cpm-normalized crosslinks from strand-agnostic 5' ends

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate deeptools

for bam in ~/clipseq/20260421/CTL/output/*.sort.dedup.bam; do
    base=$(basename "$bam" .genome.Aligned.sort.dedup.bam)
    bamCoverage --normalizeUsing CPM --Offset 1 -bs 1 -p 48 -b "$bam" -o "${base}.crosslinks.bw"
    done

for bam in ~/clipseq/20260421/OSM/output/*.sort.dedup.bam; do
    base=$(basename "$bam" .genome.Aligned.sort.dedup.bam)
    bamCoverage --normalizeUsing CPM --Offset 1 -bs 1 -p 48 -b "$bam" -o "${base}.crosslinks.bw"
    done

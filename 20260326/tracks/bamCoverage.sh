#!/bin/bash
#SBATCH --partition=cluster_long
#SBATCH --job-name=bamCoverage.sh
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

BASE_DIR=~/clipseq/20260326
OUT_DIR="${BASE_DIR}/tracks"
THREADS=52
# note this is for our clip-seq data.
# read 1 is the seq of actual RNA captured by RBP so be guided by that direction.
# when specifiying --filterRNAstrand

for condition in CTL OSM; do
    for bam in "$BASE_DIR/${condition}/output"/*.sort.dedup.bam; do
        base=$(basename "$bam" .genome.Aligned.sort.dedup.bam)

        echo "Processing: $base"

        bamCoverage \
            -b "$bam" \
            -o "${OUT_DIR}/${base}.crosslinks.plus.bw" \
            --filterRNAstrand reverse \
            --normalizeUsing CPM \
            --binSize 1 \
            --Offset 1 \
            -p "$THREADS"

        bamCoverage \
            -b "$bam" \
            -o "${OUT_DIR}/${base}.crosslinks.minus.bw" \
            --filterRNAstrand forward \
            --normalizeUsing CPM \
            --binSize 1 \
            --Offset 1 \
            -p "$THREADS"

        echo "Done: $base"
    done
done

echo "All bigWigs complete."

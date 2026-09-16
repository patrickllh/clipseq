#!/bin/bash
#SBATCH --partition=cluster_long
#SBATCH --job-name=multiBigwigSummary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

# get defined mrna regions, put together exons only, get average of crosslink signals (binsize=1)
# main output to use will be the tab file
# note: tab file coords are not typical bed; it notes down ALL single/multi-exon start/s and end/s per region per transcript in each row

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
conda activate deeptools

crosslinks_ctl=($(ls ~/clipseq/20260421/metagene_crosslinks/*Con*.crosslinks.bw | sort))
for region in 5utr cds 3utr; do
    bed=~/clipseq/20260421/extract-transcript-regions/gencode_v38_annotation_gtf_${region}.mpt.bed
    multiBigwigSummary BED-file --metagene -bs 1 -p 48 \
        --BED "$bed" \
        -b "${crosslinks_ctl[@]}" \
        -o "${region}.CTL.npz" \
        --outRawCounts "${region}.CTL.tab"
done

crosslinks_osm=($(ls ~/clipseq/20260421/metagene_crosslinks/*Sorbitol*.crosslinks.bw | sort))
for region in 5utr cds 3utr; do
    bed=~/clipseq/20260421/extract-transcript-regions/gencode_v38_annotation_gtf_${region}.mpt.bed
    multiBigwigSummary BED-file --metagene -bs 1 -p 48 \
        --BED "$bed" \
        -b "${crosslinks_osm[@]}" \
        -o "${region}.OSM.npz" \
        --outRawCounts "${region}.OSM.tab"
done

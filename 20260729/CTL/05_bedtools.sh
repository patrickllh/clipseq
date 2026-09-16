#!/bin/bash

# conda activate merge_peaks

set -euo pipefail

# ---- paths ----
TRANSCRIPT_DIR="/work/pa-hilario/clipseq/20260421/extract-transcript-regions/gtf_extract_transcript_data"
OUTPUT_DIR="/work/pa-hilario/clipseq/20260729/CTL/output"

TRANSCRIPT_BED="${TRANSCRIPT_DIR}/transcript_regions.bed"
TRANSCRIPT_SORTED="${OUTPUT_DIR}/transcript_regions.sorted.bed"

# ---- 1. Prep the transcript BED (one-time step) ----
sort -k1,1 -k2,2n "$TRANSCRIPT_BED" > "$TRANSCRIPT_SORTED"


# ---- 2. Annotate both reps ----
for rep in SNRNP70-Con-rep1_L1_R1 SNRNP70-Con-rep2_L2_R1; do

    infile="${OUTPUT_DIR}/${rep}.peaks.normed.compressed.log2FCpos.bed"
    sorted="${OUTPUT_DIR}/${rep}.peaks.sorted.bed"
    outfile="${OUTPUT_DIR}/${rep}.peaks.annotated.bed"

    sort -k1,1 -k2,2n "$infile" > "$sorted"

    bedtools intersect -a "$sorted" -b "$TRANSCRIPT_SORTED" -wa -wb -s \
        > "${OUTPUT_DIR}/${rep}.intersect.tmp"

    awk -F'\t' '{
        split($10, a, ";");
        print $1, $2, $3, $4, $5, $6, a[1], a[2], a[3]
    }' OFS='\t' "${OUTPUT_DIR}/${rep}.intersect.tmp" > "$outfile"

    rm "${OUTPUT_DIR}/${rep}.intersect.tmp"

    echo "$rep -> $outfile ($(wc -l < "$outfile") annotated peaks)"
done

#!/bin/bash
BED_DIR="$HOME/clipseq/20260421/extract-transcript-regions"

for region in 5utr cds 3utr; do
    awk 'BEGIN{OFS="\t"} {
        split($11, sizes, ",")
        total = 0
        for (i = 1; i <= $10; i++) total += sizes[i]
        print $4, total
    }' "${BED_DIR}/gencode_v38_annotation_gtf_${region}.mpt.bed" > "${region}_exonic_lengths.tsv"
done

#!/bin/bash
set -euo pipefail

OUTPUT_DIR="/work/pa-hilario/clipseq/20260729/CTL/output"

REP1="${OUTPUT_DIR}/SNRNP70-Con-rep1_L1_R1.peaks.annotated.bed"
REP2="${OUTPUT_DIR}/SNRNP70-Con-rep2_L2_R1.peaks.annotated.bed"
OUTFILE="${OUTPUT_DIR}/gene_list_union_CTL.txt"

cat \
    <(cut -f9 "$REP1" | sort -u) \
    <(cut -f9 "$REP2" | sort -u) \
    | sort -u \
    > "$OUTFILE"

n_genes=$(wc -l < "$OUTFILE")
echo "Wrote $n_genes union genes to $OUTFILE"

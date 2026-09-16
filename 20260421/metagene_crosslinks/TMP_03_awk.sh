#!/bin/bash
set -euo pipefail

# this script subtracts mV signal from SmV per respective replicate per condition
# additionally i add the transcript name based on my reference bed file

BED_DIR=~/clipseq/20260421/extract-transcript-regions # my list of transcript name and also coordinates in bed6
TAB_DIR=~/clipseq/20260421/metagene_crosslinks # crosslink cpm averaged to 1 bin size per region (excluding introns!), also with coordinates

for region in 5utr cds 3utr; do

    BED="${BED_DIR}/gencode_v38_annotation_gtf_${region}.mpt.bed"
    TAB="${TAB_DIR}/${region}.rawcounts.tab"

    # tab file structure, from deepTools multiBigwigSummary
    # Col  1=chr
    # Col  2=start
    # Col  3=end
    # Col  4=SNRNP70-Con-rep1
    # Col  5=SNRNP70-Con-rep2
    # Col  6=SNRNP70-Sorbitol-rep1
    # Col  7=SNRNP70-Sorbitol-rep2
    # Col  8=mVenus-Con-rep5
    # Col  9=mVenus-Con-rep6
    # Col  10=mVenus-Sorbitol-rep5
    # Col  11=mVenus-Sorbitol-rep6
    # then i join based on concatenated chr_start_end to ensure unique rows

     join -t $'\t' -1 1 -2 1 \
        <(awk 'BEGIN{OFS="\t"} {print $1"_"$2"_"$3, $1, $2, $3, $4}' "$BED" | sort -k1,1) \
        <(tail -n +2 "$TAB" | awk 'BEGIN{OFS="\t"} {print $1"_"$2"_"$3, $4-$8, $5-$9}' | sort -k1,1) \
     | awk 'BEGIN{OFS="\t"; print "transcript_id\tchr\tstart\tend\tCTL_rep1_SNRNP70_minus_mVenus\tCTL_rep2_SNRNP70_minus_mVenus"} \
           {print $2, $3, $4, $5, $6, $7}' \
        > "${region}.CTL.tsv"

     join -t $'\t' -1 1 -2 1 \
        <(awk 'BEGIN{OFS="\t"} {print $1"_"$2"_"$3, $1, $2, $3, $4}' "$BED" | sort -k1,1) \
        <(tail -n +2 "$TAB" | awk 'BEGIN{OFS="\t"} {print $1"_"$2"_"$3, $6-$10, $7-$11}' | sort -k1,1) \
     | awk 'BEGIN{OFS="\t"; print "transcript_id\tchr\tstart\tend\tOSM_rep1_SNRNP70_minus_mVenus\tOSM_rep2_SNRNP70_minus_mVenus"} \
           {print $2, $3, $4, $5, $6, $7}' \
        > "${region}.OSM.tsv"
done

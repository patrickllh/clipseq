#!/bin/bash

# Please set-up first!

# Set-up tool for extracting 1 transcript per gene
conda create -n rbpbench -c conda-forge -c bioconda rbpbench
conda activate rbpbench

# Set-up tool for extracting mRNA regions
cd ~
git clone https://github.com/stephenfloor/extract-transcript-regions.git
cd extract-transcript-regions
chmod +x extract_transcript_regions.py

# Set-up GRCh38 genome & GENCODE v38 annotation
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz && \
gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz && \
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz && \
gunzip gencode.v38.annotation.gtf.gz

# mRNA annotation starts here:
# RBPBench: extract 1 transcript per gene as BED12, based on "basic" tag, TSL, length; skip TEC transcripts!
gtf_extract_transcript_data.py \
--gtf gencode.v38.annotation.gtf \
--genome GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--out gtf_extract_transcript_data \
--mrna-only \
--skip-tec

# extract-transcript-regions: extract 5UTR, CDS, 3UTR, introns, etc., per transcript as BED12
python extract_transcript_regions.py -i gencode.v38.annotation.gtf -o gencode_v38_annotation_gtf --gtf

# filter mRNA regions into only transcripts selected by RBPBench
for f in \
gencode_v38_annotation_gtf_5utr.bed \
gencode_v38_annotation_gtf_cds.bed \
gencode_v38_annotation_gtf_3utr.bed
do
awk -F'\t' '
NR==FNR{
    split($4,a,";");
    keep[a[1]]=1;
    next
}
{
    id=$4
    sub(/_(5utr|3utr|cds)$/, "", id)
    if (id in keep) print
}' gtf_extract_transcript_data/transcript_regions.bed "$f" > "${f%.bed}.mpt.bed"
done

# Check: union of transcripts from mRNA regions `.mpt.bed` was the same as those from `transcript_regions.bed`.

# get list of transcripts in transcript_regions.bed
cut -f4 gtf_extract_transcript_data/transcript_regions.bed | cut -d';' -f1 | sort -u > ref.ids

# get union of transcripts from the *mpt.bed mRNA regions
cut -f4 *.mpt.bed \
| sed -E 's/_(5utr|3utr|cds)$//' \
| sort -u > mpt.ids

# identify transcripts found only in transcript_regions.bed
comm -23 ref.ids mpt.ids # empty

# identify transcripts found only in *mpt.bed
comm -13 ref.ids mpt.ids | head # empty

# filter intron regions as BED6 file into only transcripts selected by RBPBench
awk '
BEGIN{FS=OFS="\t"}
FNR==NR {
    split($4,a,";");
    keep[a[1]]=1;
    next
}
{
    split($4,b,";");
    tx=b[1];
    feature=b[4];

    if (feature=="intron" && tx in keep)
        print
}
' gtf_extract_transcript_data/transcript_regions.bed gtf_extract_transcript_data/exon_intron_regions.bed > intron_regions.mpt.bed

# quick additional sanity check to see if feature-containing transcripts are below total transcripts
wc -l gencode_v38_annotation_gtf_*.mpt.bed
#   19918 gencode_v38_annotation_gtf_3utr.mpt.bed
#   19505 gencode_v38_annotation_gtf_5utr.mpt.bed
#   20413 gencode_v38_annotation_gtf_cds.mpt.bed
cut -f4 intron_regions.mpt.bed | sort | uniq | wc -l
#   19123
wc -l gtf_extract_transcript_data/transcript_regions.bed
#   20413 gtf_extract_transcript_data/transcript_regions.bed
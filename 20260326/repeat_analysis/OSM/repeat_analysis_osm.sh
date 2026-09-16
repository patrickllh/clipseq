#!/bin/bash

# download RepeatMasker for hg38
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz ./repeatmasker.grch38.tsv.gz

# RepeatMasker data structure
zcat repeatmasker.grch38.tsv.gz | head -n 2
# 585     463     13      6       17      chr1    10000   10468   -248945954      +       (TAACCC)n     Simple_repeat    Simple_repeat   1       471     0       1
# 585     3612    114     215     13      chr1    10468   11447   -248944975      -       TAR1    Satellite      telo    -399    1712    483     2

# regions per repeat family, which we can survey
zcat repeatmasker.grch38.tsv.gz | cut -f12 | sort | uniq -c | sort -nr
# 1910631 SINE
# 1614481 LINE
#  770551 LTR
#  724562 Simple_repeat
#  512404 DNA
#  106053 Low_complexity
#    9133 Satellite
#    5974 Retroposon
#    5777 LTR?
#    5732 Unknown
#    4686 snRNA
#    3364 DNA?
#    2164 tRNA
#    1953 rRNA
#    1820 RC
#    1745 srpRNA
#    1484 scRNA
#     721 RNA
#     417 RC?
#      38 SINE?

# convert RepeatMasker data to BED6
zcat ~/skipper/annotations/repeatmasker.grch38.tsv.gz \
| awk 'BEGIN{OFS="\t"} {print $6,$7,$8,$11,$12,$10}' \
> repeatmasker.grch38.bed

# how many genomic intervals for repeats do we get? 
zcat repeatmasker.grch38.tsv.gz | wc -l # 5683690 repeats
wc -l repeatmasker.grch38.bed # 5683690 repeats, just sanity check


# we now inspect mapping profile of a STAR mapping run which allowed multimapping reads (each run is different/not deterministic)
for log in ~/clipseq/20260326/OSM/output/*.Log.final.out; do echo $log; cat $log; done

# /work/pa-hilario/clipseq/20260326/OSM/output/mVenus-Sorbitol-rep5_L1_R1.genome.Log.final.out
#
#                           Number of input reads |       925360
#                       Average input read length |       62
#
#                                     UNIQUE READS:
#                    Uniquely mapped reads number |       112054
#                         Uniquely mapped reads % |       12.11%
#                           Average mapped length |       52.33
#
#
#                              MULTI-MAPPING READS:
#         Number of reads mapped to multiple loci |       350900
#              % of reads mapped to multiple loci |       37.92%
#         Number of reads mapped to too many loci |       20
#              % of reads mapped to too many loci |       0.00%
#
#                                   UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |       13704
#        % of reads unmapped: too many mismatches |       1.48%
#             Number of reads unmapped: too short |       383063
#                  % of reads unmapped: too short |       41.40%
#                 Number of reads unmapped: other |       65619
#                      % of reads unmapped: other |       7.09%
#
# /work/pa-hilario/clipseq/20260326/OSM/output/mVenus-Sorbitol-rep6_L2_R1.genome.Log.final.out
#
#                           Number of input reads |       2400846
#                       Average input read length |       87
#
#                                     UNIQUE READS:
#                    Uniquely mapped reads number |       238439
#                         Uniquely mapped reads % |       9.93%
#                           Average mapped length |       82.38
#
#                              MULTI-MAPPING READS:
#         Number of reads mapped to multiple loci |       1792429
#              % of reads mapped to multiple loci |       74.66%
#         Number of reads mapped to too many loci |       70
#              % of reads mapped to too many loci |       0.00%
#
#                                   UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |       56764
#        % of reads unmapped: too many mismatches |       2.36%
#             Number of reads unmapped: too short |       300417
#                  % of reads unmapped: too short |       12.51%
#                 Number of reads unmapped: other |       12727
#                      % of reads unmapped: other |       0.53%
#
# /work/pa-hilario/clipseq/20260326/OSM/output/SNRNP70-Sorbitol-rep1_L1_R1.genome.Log.final.out
#
#                           Number of input reads |       9144479
#                       Average input read length |       62
#
#                                     UNIQUE READS:
#                    Uniquely mapped reads number |       514712
#                         Uniquely mapped reads % |       5.63%
#                           Average mapped length |       61.46
#
#                              MULTI-MAPPING READS:
#         Number of reads mapped to multiple loci |       8040206
#              % of reads mapped to multiple loci |       87.92%
#         Number of reads mapped to too many loci |       93
#              % of reads mapped to too many loci |       0.00%
#
#                                   UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |       71059
#        % of reads unmapped: too many mismatches |       0.78%
#             Number of reads unmapped: too short |       506073
#                  % of reads unmapped: too short |       5.53%
#                 Number of reads unmapped: other |       12336
#                      % of reads unmapped: other |       0.13%
#
# /work/pa-hilario/clipseq/20260326/OSM/output/SNRNP70-Sorbitol-rep2_L2_R1.genome.Log.final.out
#
#                           Number of input reads |       9216267
#                       Average input read length |       65
#
#                                     UNIQUE READS:
#                    Uniquely mapped reads number |       544886
#                         Uniquely mapped reads % |       5.91%
#                           Average mapped length |       64.22
#
#                              MULTI-MAPPING READS:
#         Number of reads mapped to multiple loci |       7960900
#              % of reads mapped to multiple loci |       86.38%
#         Number of reads mapped to too many loci |       251
#              % of reads mapped to too many loci |       0.00%
#
#                                   UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |       75759
#        % of reads unmapped: too many mismatches |       0.82%
#             Number of reads unmapped: too short |       624021
#                  % of reads unmapped: too short |       6.77%
#                 Number of reads unmapped: other |       10450
#                      % of reads unmapped: other |       0.11%

# most SNRNP70 reads were multimapping (rep1: 87.92%, rep2: 86.38%)
# some uniquely mapped SNRNP70 reads exist (rep1: 5.63%, rep2: 5.91%) with read length ~61-64nt
# note: that these are not yet deduplicated. from here-on I use deduplicated reads, which don't have a report

# to check deduplication rate, I compare total mapped reads pre- and post-deduplication

# pre-deduplication
for bam in ~/clipseq/20260326/OSM/output/*.genome.Aligned.out.bam; do echo $bam; samtools flagstat $bam; done
# /work/pa-hilario/clipseq/20260326/OSM/output/mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.out.bam
# 925360 + 0 in total (QC-passed reads + QC-failed reads)
# 925360 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 462954 + 0 mapped (50.03% : N/A)
# 462954 + 0 primary mapped (50.03% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# /work/pa-hilario/clipseq/20260326/OSM/output/mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.out.bam
# 2400846 + 0 in total (QC-passed reads + QC-failed reads)
# 2400846 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 2030868 + 0 mapped (84.59% : N/A)
# 2030868 + 0 primary mapped (84.59% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# /work/pa-hilario/clipseq/20260326/OSM/output/SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.out.bam
# 9144479 + 0 in total (QC-passed reads + QC-failed reads)
# 9144479 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 8554918 + 0 mapped (93.55% : N/A)
# 8554918 + 0 primary mapped (93.55% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# /work/pa-hilario/clipseq/20260326/OSM/output/SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.out.bam
# 9216267 + 0 in total (QC-passed reads + QC-failed reads)
# 9216267 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 8505786 + 0 mapped (92.29% : N/A)
# 8505786 + 0 primary mapped (92.29% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# post-deduplication
for bam in ~/clipseq/20260326/OSM/output/*.genome.Aligned.sort.dedup.bam; do echo $bam; samtools flagstat $bam; done

# /work/pa-hilario/clipseq/20260326/OSM/output/mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.sort.dedup.bam
# 153905 + 0 in total (QC-passed reads + QC-failed reads)
# 153905 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 153905 + 0 mapped (100.00% : N/A)
# 153905 + 0 primary mapped (100.00% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
#
# /work/pa-hilario/clipseq/20260326/OSM/output/mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.sort.dedup.bam
# 1382898 + 0 in total (QC-passed reads + QC-failed reads)
# 1382898 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 1382898 + 0 mapped (100.00% : N/A)
# 1382898 + 0 primary mapped (100.00% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
#
# /work/pa-hilario/clipseq/20260326/OSM/output/SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.sort.dedup.bam
# 2492229 + 0 in total (QC-passed reads + QC-failed reads)
# 2492229 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 2492229 + 0 mapped (100.00% : N/A)
# 2492229 + 0 primary mapped (100.00% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
#
# /work/pa-hilario/clipseq/20260326/OSM/output/SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.sort.dedup.bam
# 1172483 + 0 in total (QC-passed reads + QC-failed reads)
# 1172483 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 1172483 + 0 mapped (100.00% : N/A)
# 1172483 + 0 primary mapped (100.00% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# the deduplication rate resultss are as follows
# | Sample                | Pre-dedup mapped | Post-dedup mapped | Retained % | Deduplicated % |
# | --------------------- | ---------------: | ----------------: | ---------: | -------------: |
# | mVenus-Sorbitol-rep5  |          462,954 |           153,905 |     33.24% |     **66.76%** |
# | mVenus-Sorbitol-rep6  |        2,030,868 |         1,382,898 |     68.09% |     **31.91%** |
# | SNRNP70-Sorbitol-rep1 |        8,554,918 |         2,492,229 |     29.13% |     **70.87%** |
# | SNRNP70-Sorbitol-rep2 |        8,505,786 |         1,172,483 |     13.78% |     **86.22%** |

# deduplication was quite drastic (70%-86% for SNRNP70), leaving only 1-2 million reads
# if library complexity is low (e.g. due to U1 snRNA) then perhaps this isn't surprising

# to survey multimapping reads, I filter bams for multimapping reads only, via "NH:i:" >1 in their read headers
conda activate skipper

# SNRNP70 rep1 multimappers
samtools view -h ~/clipseq/20260326/OSM/output/SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.sort.dedup.bam \
| awk 'BEGIN{OFS="\t"}
/^@/ {print; next}
{
  nh=1
  for(i=12;i<=NF;i++){
    if($i ~ /^NH:i:/){
      split($i,a,":")
      nh=a[3]
    }
  }
  if(nh>1) print
}' \
| samtools view -b -o SNRNP70_OSM_1.multimap.bam

# SNRNP70 rep2 multimappers
samtools view -h ~/clipseq/20260326/OSM/output/SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.sort.dedup.bam \
| awk 'BEGIN{OFS="\t"}
/^@/ {print; next}
{
  nh=1
  for(i=12;i<=NF;i++){
    if($i ~ /^NH:i:/){
      split($i,a,":")
      nh=a[3]
    }
  }
  if(nh>1) print
}' \
| samtools view -b -o SNRNP70_OSM_2.multimap.bam

# mVenus rep1 multimappers
samtools view -h ~/clipseq/20260326/OSM/output/mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.sort.dedup.bam \
| awk 'BEGIN{OFS="\t"}
/^@/ {print; next}
{
  nh=1
  for(i=12;i<=NF;i++){
    if($i ~ /^NH:i:/){
      split($i,a,":")
      nh=a[3]
    }
  }
  if(nh>1) print
}' \
| samtools view -b -o mVenus_OSM_1.multimap.bam

# mVenus rep2 multimappers
samtools view -h ~/clipseq/20260326/OSM/output/mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.sort.dedup.bam \
| awk 'BEGIN{OFS="\t"}
/^@/ {print; next}
{
  nh=1
  for(i=12;i<=NF;i++){
    if($i ~ /^NH:i:/){
      split($i,a,":")
      nh=a[3]
    }
  }
  if(nh>1) print
}' \
| samtools view -b -o mVenus_OSM_2.multimap.bam

# find reads which have at least 70% overlaps to repeat list

# SNRNP70 rep1
bedtools intersect \
-a SNRNP70_OSM_1.multimap.bam \
-b repeatmasker.grch38.bed  \
-bed -wa -wb -s -f 0.7 > SNRNP70_OSM_1.multimap.inrepeats.tsv

# SNRNP70 rep2
bedtools intersect \
-a SNRNP70_OSM_2.multimap.bam \
-b repeatmasker.grch38.bed  \
-bed -wa -wb -s -f 0.7 > SNRNP70_OSM_2.multimap.inrepeats.tsv

# mVenus rep1
bedtools intersect \
-a mVenus_OSM_1.multimap.bam \
-b repeatmasker.grch38.bed  \
-bed -wa -wb -s -f 0.7 > mVenus_OSM_1.multimap.inrepeats.tsv

# mVenus rep2
bedtools intersect \
-a mVenus_OSM_2.multimap.bam \
-b repeatmasker.grch38.bed  \
-bed -wa -wb -s -f 0.7 > mVenus_OSM_2.multimap.inrepeats.tsv

# number of multimapping reads (deduplicated)
for multimapper in *.multimap.bam; do echo $multimapper; samtools view $multimapper | wc -l; done
# mVenus_OSM_1.multimap.bam
# 127952
# mVenus_OSM_2.multimap.bam
# 1278289
# SNRNP70_OSM_1.multimap.bam
# 2092180
# SNRNP70_OSM_2.multimap.bam
# 1101163

# number of multimapping reads in repeats
wc -l *.multimap.inrepeats.tsv
#    107706 mVenus_OSM_1.multimap.inrepeats.tsv
#   1165506 mVenus_OSM_2.multimap.inrepeats.tsv
#   1928611 SNRNP70_OSM_1.multimap.inrepeats.tsv
#   1040749 SNRNP70_OSM_2.multimap.inrepeats.tsv

# SNRNP70 multimapping reads mostly map to repeats (rep1: 92%, rep2: 94%)
# mVenus multimapping reads also mostly map to repeats (rep1: 85%, rep2: 91%)

# importantly, we check the overall repeat-binding profile of multimappers per library
# note: sometimes 1 read can map to >1 repeat region if parts of the read intersect with multiple regions

for tsv in *.multimap.inrepeats.tsv; do echo $tsv; cut -f17 $tsv | sort | uniq -c | sort -nr
# mVenus_OSM_1.multimap.inrepeats.tsv
#   95991 rRNA
#    3894 snRNA
#    3749 tRNA
#    1167 LINE
#     917 SINE
#     872 Simple_repeat
#     544 LTR
#     207 DNA
#     123 srpRNA
#     103 Low_complexity
#      46 Satellite
#      39 Retroposon
#      29 scRNA
#      25 RNA

# mVenus_OSM_2.multimap.inrepeats.tsv
# 1140464 rRNA
#    7793 tRNA
#    7554 snRNA
#    3087 Simple_repeat
#    1498 DNA
#    1481 LINE
#    1357 SINE
#     863 srpRNA
#     828 LTR
#     227 RNA
#     182 Low_complexity
#      74 Satellite
#      60 scRNA
#      35 Retroposon
#       2 Unknown
#       1 LTR?

# SNRNP70_OSM_1.multimap.inrepeats.tsv
#  986044 snRNA
#  906132 rRNA
#   13135 tRNA
#    9645 Simple_repeat
#    4297 DNA
#    3893 RNA
#    2723 SINE
#     851 LTR
#     655 LINE
#     415 Low_complexity
#     376 Retroposon
#     304 srpRNA
#      99 Satellite
#      36 scRNA
#       3 LTR?
#       2 RC
#       1 Unknown

# SNRNP70_OSM_2.multimap.inrepeats.tsv
#  603816 snRNA
#  420462 rRNA
#    6084 tRNA
#    4117 Simple_repeat
#    2302 DNA
#    1250 SINE
#     933 RNA
#     547 LTR
#     518 LINE
#     387 Retroposon
#     155 Low_complexity
#     102 srpRNA
#      58 Satellite
#      17 scRNA
#       1 LTR?

# how many U1 snRNA-mapping multimapping reads?
for tsv in *.multimap.inrepeats.tsv; do echo $tsv; grep snRNA $tsv | cut -f16 | grep -w "U1" | wc -l; done
# mVenus_OSM_1.multimap.inrepeats.tsv
# 2411
# mVenus_OSM_2.multimap.inrepeats.tsv
# 5155
# SNRNP70_OSM_1.multimap.inrepeats.tsv
# 984683
# SNRNP70_OSM_2.multimap.inrepeats.tsv
# 603152

# SNRNP70 OSM libraries mostly have multimapping reads, with around half mapping to snRNAs (rep1: 47% rep2: 54%), and another to rRNAs (rep1: 43% rep2: 38%)
# importantly, SNRNP70 multimapping snRNA-mapped reads were predominated with U1 snRNA (rep1: 99% rep2: 99%)
# in other words, half of them in SNRNP70 are virtually U1 snRNA, while half are predominantly rRNAs among other repeat sequences.
# moreover, repeating from above, these multimapping reads are mostly found in repeats (rep1: 92%, rep2: 94%)
# taken together, SNRNP70 mostly binds to repeat regions, half of which are U1 snRNAs, and the other half predominantly rRNAs.

# in contrast, mVenus multimappers were mostly rRNAs (rep1: 75%, rep2: 89%).
# these were followed either by tRNAs (rep1: 3%, rep2: 0.6%) or snRNAs (rep1: 3%, rep2: 0.5%)

# note again, that there is a possibility of reads mapping to multiple regions if part of them maps to one and the other
# however it seems to be just for a very small number of reads, relatively:
for tsv in *.multimap.inrepeats.tsv; do echo $tsv; cut -f4 $tsv | sort | uniq -d | wc -l; done
# mVenus_OSM_1.multimap.inrepeats.tsv
# 6 (out of 2411)
# mVenus_OSM_2.multimap.inrepeats.tsv
# 4 (out of 5155)
# SNRNP70_OSM_1.multimap.inrepeats.tsv
# 10 (out of 984683)
# SNRNP70_OSM_2.multimap.inrepeats.tsv
# 10 (out of 603152)
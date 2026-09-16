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
zcat repeatmasker.grch38.tsv.gz \
| awk 'BEGIN{OFS="\t"} {print $6,$7,$8,$11,$12,$10}' \
> repeatmasker.grch38.bed

# how many genomic intervals for repeats do we get? 
zcat repeatmasker.grch38.tsv.gz | wc -l # 5683690 repeats
wc -l repeatmasker.grch38.bed # 5683690 repeats, just sanity check


# we now inspect mapping profile of a STAR mapping run which allowed multimapping reads (each run is different/not deterministic)
for log in ~/clipseq/20260326/CTL/output/*.Log.final.out; do echo $log; cat $log; done
# /work/pa-hilario/clipseq/20260326/CTL/output/mVenus-Con-rep5_L1_R1.genome.Log.final.out
#                                  Started job on |       Mar 27 15:26:47
#                              Started mapping on |       Mar 27 15:27:18
#                                     Finished on |       Mar 27 15:27:41
#        Mapping speed, Million of reads per hour |       79.00

#                           Number of input reads |       504732
#                       Average input read length |       58
#                                     UNIQUE READS:
#                    Uniquely mapped reads number |       50099
#                         Uniquely mapped reads % |       9.93%
#                           Average mapped length |       43.54
#                        Number of splices: Total |       2519
#             Number of splices: Annotated (sjdb) |       0
#                        Number of splices: GT/AG |       2021
#                        Number of splices: GC/AG |       119
#                        Number of splices: AT/AC |       5
#                Number of splices: Non-canonical |       374
#                       Mismatch rate per base, % |       3.13%
#                          Deletion rate per base |       0.06%
#                         Deletion average length |       2.50
#                         Insertion rate per base |       0.03%
#                        Insertion average length |       1.24
#                              MULTI-MAPPING READS:
#         Number of reads mapped to multiple loci |       365782
#              % of reads mapped to multiple loci |       72.47%
#         Number of reads mapped to too many loci |       155
#              % of reads mapped to too many loci |       0.03%
#                                   UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |       7343
#        % of reads unmapped: too many mismatches |       1.45%
#             Number of reads unmapped: too short |       77084
#                  % of reads unmapped: too short |       15.27%
#                 Number of reads unmapped: other |       4269
#                      % of reads unmapped: other |       0.85%
#                                   CHIMERIC READS:
#                        Number of chimeric reads |       0
#                             % of chimeric reads |       0.00%
# /work/pa-hilario/clipseq/20260326/CTL/output/mVenus-Con-rep6_L2_R1.genome.Log.final.out
#                                  Started job on |       Mar 27 15:26:47
#                              Started mapping on |       Mar 27 15:27:18
#                                     Finished on |       Mar 27 15:27:46
#        Mapping speed, Million of reads per hour |       250.17

#                           Number of input reads |       1945761
#                       Average input read length |       87
#                                     UNIQUE READS:
#                    Uniquely mapped reads number |       187664
#                         Uniquely mapped reads % |       9.64%
#                           Average mapped length |       73.74
#                        Number of splices: Total |       15589
#             Number of splices: Annotated (sjdb) |       0
#                        Number of splices: GT/AG |       11906
#                        Number of splices: GC/AG |       305
#                        Number of splices: AT/AC |       24
#                Number of splices: Non-canonical |       3354
#                       Mismatch rate per base, % |       1.59%
#                          Deletion rate per base |       0.05%
#                         Deletion average length |       1.70
#                         Insertion rate per base |       0.01%
#                        Insertion average length |       1.10
#                              MULTI-MAPPING READS:
#         Number of reads mapped to multiple loci |       1343279
#              % of reads mapped to multiple loci |       69.04%
#         Number of reads mapped to too many loci |       88
#              % of reads mapped to too many loci |       0.00%
#                                   UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |       65170
#        % of reads unmapped: too many mismatches |       3.35%
#             Number of reads unmapped: too short |       327798
#                  % of reads unmapped: too short |       16.85%
#                 Number of reads unmapped: other |       21762
#                      % of reads unmapped: other |       1.12%
#                                   CHIMERIC READS:
#                        Number of chimeric reads |       0
#                             % of chimeric reads |       0.00%
# /work/pa-hilario/clipseq/20260326/CTL/output/SNRNP70-Con-rep1_L1_R1.genome.Log.final.out
#                                  Started job on |       Mar 27 15:26:47
#                              Started mapping on |       Mar 27 15:27:19
#                                     Finished on |       Mar 27 15:28:00
#        Mapping speed, Million of reads per hour |       720.36

#                           Number of input reads |       8204139
#                       Average input read length |       63
#                                     UNIQUE READS:
#                    Uniquely mapped reads number |       345115
#                         Uniquely mapped reads % |       4.21%
#                           Average mapped length |       56.63
#                        Number of splices: Total |       20294
#             Number of splices: Annotated (sjdb) |       0
#                        Number of splices: GT/AG |       16518
#                        Number of splices: GC/AG |       385
#                        Number of splices: AT/AC |       138
#                Number of splices: Non-canonical |       3253
#                       Mismatch rate per base, % |       2.08%
#                          Deletion rate per base |       0.20%
#                         Deletion average length |       2.62
#                         Insertion rate per base |       0.01%
#                        Insertion average length |       1.75
#                              MULTI-MAPPING READS:
#         Number of reads mapped to multiple loci |       7505407
#              % of reads mapped to multiple loci |       91.48%
#         Number of reads mapped to too many loci |       119
#              % of reads mapped to too many loci |       0.00%
#                                   UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |       45596
#        % of reads unmapped: too many mismatches |       0.56%
#             Number of reads unmapped: too short |       299134
#                  % of reads unmapped: too short |       3.65%
#                 Number of reads unmapped: other |       8768
#                      % of reads unmapped: other |       0.11%
#                                   CHIMERIC READS:
#                        Number of chimeric reads |       0
#                             % of chimeric reads |       0.00%
# /work/pa-hilario/clipseq/20260326/CTL/output/SNRNP70-Con-rep2_L2_R1.genome.Log.final.out
#                                  Started job on |       Mar 27 15:26:47
#                              Started mapping on |       Mar 27 15:27:18
#                                     Finished on |       Mar 27 15:27:58
#        Mapping speed, Million of reads per hour |       695.21

#                           Number of input reads |       7724609
#                       Average input read length |       62
#                                     UNIQUE READS:
#                    Uniquely mapped reads number |       286999
#                         Uniquely mapped reads % |       3.72%
#                           Average mapped length |       55.87
#                        Number of splices: Total |       13113
#             Number of splices: Annotated (sjdb) |       0
#                        Number of splices: GT/AG |       12511
#                        Number of splices: GC/AG |       269
#                        Number of splices: AT/AC |       32
#                Number of splices: Non-canonical |       301
#                       Mismatch rate per base, % |       1.61%
#                          Deletion rate per base |       0.09%
#                         Deletion average length |       3.02
#                         Insertion rate per base |       0.01%
#                        Insertion average length |       1.63
#                              MULTI-MAPPING READS:
#         Number of reads mapped to multiple loci |       7189560
#              % of reads mapped to multiple loci |       93.07%
#         Number of reads mapped to too many loci |       78
#              % of reads mapped to too many loci |       0.00%
#                                   UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |       32650
#        % of reads unmapped: too many mismatches |       0.42%
#             Number of reads unmapped: too short |       207018
#                  % of reads unmapped: too short |       2.68%
#                 Number of reads unmapped: other |       8304
#                      % of reads unmapped: other |       0.11%
#                                   CHIMERIC READS:
#                        Number of chimeric reads |       0
#                             % of chimeric reads |       0.00%
# most SNRNP70 reads were multimapping (rep1: 91.48%, rep2: 93.07%)
# some uniquely mapped SNRNP70 reads exist (rep1: 4.21%, rep2: 3.72%) with read length ~62-63nt
# note: that these are not yet deduplicated. from here-on I use deduplicated reads, which don't have a report

# to check deduplication rate, I compare total mapped reads pre- and post-deduplication

# pre-deduplication
for bam in ~/clipseq/20260326/CTL/output/*.genome.Aligned.out.bam; do echo $bam; samtools flagstat $bam; done
# /work/pa-hilario/clipseq/20260326/CTL/output/mVenus-Con-rep5_L1_R1.genome.Aligned.out.bam
# 504732 + 0 in total (QC-passed reads + QC-failed reads)
# 504732 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 415881 + 0 mapped (82.40% : N/A)
# 415881 + 0 primary mapped (82.40% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
# /work/pa-hilario/clipseq/20260326/CTL/output/mVenus-Con-rep6_L2_R1.genome.Aligned.out.bam
# 1945761 + 0 in total (QC-passed reads + QC-failed reads)
# 1945761 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 1530943 + 0 mapped (78.68% : N/A)
# 1530943 + 0 primary mapped (78.68% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
# /work/pa-hilario/clipseq/20260326/CTL/output/SNRNP70-Con-rep1_L1_R1.genome.Aligned.out.bam
# 8204139 + 0 in total (QC-passed reads + QC-failed reads)
# 8204139 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 7850522 + 0 mapped (95.69% : N/A)
# 7850522 + 0 primary mapped (95.69% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
# /work/pa-hilario/clipseq/20260326/CTL/output/SNRNP70-Con-rep2_L2_R1.genome.Aligned.out.bam
# 7724609 + 0 in total (QC-passed reads + QC-failed reads)
# 7724609 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 7476559 + 0 mapped (96.79% : N/A)
# 7476559 + 0 primary mapped (96.79% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# post-deduplication
for bam in ~/clipseq/20260326/CTL/output/*.genome.Aligned.sort.dedup.bam; do echo $bam; samtools flagstat $bam; done
# /work/pa-hilario/clipseq/20260326/CTL/output/mVenus-Con-rep5_L1_R1.genome.Aligned.sort.dedup.bam
# 78393 + 0 in total (QC-passed reads + QC-failed reads)
# 78393 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 78393 + 0 mapped (100.00% : N/A)
# 78393 + 0 primary mapped (100.00% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
# /work/pa-hilario/clipseq/20260326/CTL/output/mVenus-Con-rep6_L2_R1.genome.Aligned.sort.dedup.bam
# 195523 + 0 in total (QC-passed reads + QC-failed reads)
# 195523 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 195523 + 0 mapped (100.00% : N/A)
# 195523 + 0 primary mapped (100.00% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
# /work/pa-hilario/clipseq/20260326/CTL/output/SNRNP70-Con-rep1_L1_R1.genome.Aligned.sort.dedup.bam
# 1586439 + 0 in total (QC-passed reads + QC-failed reads)
# 1586439 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 1586439 + 0 mapped (100.00% : N/A)
# 1586439 + 0 primary mapped (100.00% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
# /work/pa-hilario/clipseq/20260326/CTL/output/SNRNP70-Con-rep2_L2_R1.genome.Aligned.sort.dedup.bam
# 2253025 + 0 in total (QC-passed reads + QC-failed reads)
# 2253025 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 2253025 + 0 mapped (100.00% : N/A)
# 2253025 + 0 primary mapped (100.00% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# the deduplication rate results are as follows
# | Sample               | Pre-dedup mapped | Post-dedup mapped | Retained % | Deduplicated % |
# | -------------------- | ---------------: | ----------------: | ---------: | -------------: |
# | mVenus-Con-rep5      |          415,881 |            78,393 |     18.85% |     **81.15%** |
# | mVenus-Con-rep6      |        1,530,943 |           195,523 |     12.77% |     **87.23%** |
# | SNRNP70-Con-rep1     |        7,850,522 |         1,586,439 |     20.21% |     **79.79%** |
# | SNRNP70-Con-rep2     |        7,476,559 |         2,253,025 |     30.13% |     **69.87%** |

# deduplication was quite drastic (70%-86% for SNRNP70), leaving only 1-2 million reads
# if library complexity is low (e.g. due to U1 snRNA) then perhaps this isn't surprising

# to survey multimapping reads, I filter bams for multimapping reads only, via "NH:i:" >1 in their read headers
conda activate skipper

# SNRNP70 rep1 multimappers
samtools view -h ~/clipseq/20260326/CTL/output/SNRNP70-Con-rep1_L1_R1.genome.Aligned.sort.dedup.bam \
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
| samtools view -b -o SNRNP70_CTL_1.multimap.bam

# SNRNP70 rep2 multimappers
samtools view -h ~/clipseq/20260326/CTL/output/SNRNP70-Con-rep2_L2_R1.genome.Aligned.sort.dedup.bam \
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
| samtools view -b -o SNRNP70_CTL_2.multimap.bam

# mVenus rep1 multimappers
samtools view -h ~/clipseq/20260326/CTL/output/mVenus-Con-rep5_L1_R1.genome.Aligned.sort.dedup.bam \
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
| samtools view -b -o mVenus_CTL_1.multimap.bam

# mVenus rep2 multimappers
samtools view -h ~/clipseq/20260326/CTL/output/mVenus-Con-rep6_L2_R1.genome.Aligned.sort.dedup.bam \
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
| samtools view -b -o mVenus_CTL_2.multimap.bam

# find reads which have at least 70% overlaps to repeat list

# SNRNP70 rep1
bedtools intersect \
-a SNRNP70_CTL_1.multimap.bam \
-b repeatmasker.grch38.bed  \
-bed -wa -wb -s -f 0.7 > SNRNP70_CTL_1.multimap.inrepeats.tsv

# SNRNP70 rep2
bedtools intersect \
-a SNRNP70_CTL_2.multimap.bam \
-b repeatmasker.grch38.bed  \
-bed -wa -wb -s -f 0.7 > SNRNP70_CTL_2.multimap.inrepeats.tsv

# mVenus rep1
bedtools intersect \
-a mVenus_CTL_1.multimap.bam \
-b repeatmasker.grch38.bed  \
-bed -wa -wb -s -f 0.7 > mVenus_CTL_1.multimap.inrepeats.tsv

# mVenus rep2
bedtools intersect \
-a mVenus_CTL_2.multimap.bam \
-b repeatmasker.grch38.bed  \
-bed -wa -wb -s -f 0.7 > mVenus_CTL_2.multimap.inrepeats.tsv

# number of multimapping reads (deduplicated)
for multimapper in *.multimap.bam; do echo $multimapper; samtools view $multimapper | wc -l; done
# mVenus_CTL_1.multimap.bam
# 66381
# mVenus_CTL_2.multimap.bam
# 178396
# SNRNP70_CTL_1.multimap.bam
# 1473216
# SNRNP70_CTL_2.multimap.bam
# 2035794

# number of multimapping reads in repeats
wc -l *.multimap.inrepeats.tsv
  #   53804 mVenus_CTL_1.multimap.inrepeats.tsv
  #  153626 mVenus_CTL_2.multimap.inrepeats.tsv
  # 1399508 SNRNP70_CTL_1.multimap.inrepeats.tsv
  # 1934095 SNRNP70_CTL_2.multimap.inrepeats.tsv

# SNRNP70 multimapping reads mostly map to repeats (rep1: 95%, rep2: 95%)
# mVenus multimapping reads also mostly map to repeats (rep1: 81%, rep2: 86%)

# importantly, we check the overall repeat-binding profile of multimappers per library
# note: sometimes 1 read can map to >1 repeat region if parts of the read intersect with multiple regions

for tsv in *.multimap.inrepeats.tsv; do echo $tsv; cut -f17 $tsv | sort | uniq -c | sort -nr; done
# mVenus_CTL_1.multimap.inrepeats.tsv
#   48247 rRNA
#    1399 snRNA
#    1336 tRNA
#    1130 Simple_repeat
#     624 SINE
#     456 LINE
#     311 LTR
#     168 DNA
#      65 Low_complexity
#      28 scRNA
#      20 Satellite
#       9 RNA
#       8 srpRNA
#       3 Retroposon
# mVenus_CTL_2.multimap.inrepeats.tsv
#  146568 rRNA
#    2885 snRNA
#     923 tRNA
#     898 Simple_repeat
#     708 LINE
#     537 SINE
#     462 DNA
#     432 LTR
#      99 srpRNA
#      38 Retroposon
#      27 RNA
#      25 scRNA
#      12 Satellite
#       8 Low_complexity
#       2 LTR?
#       2 DNA?
# SNRNP70_CTL_1.multimap.inrepeats.tsv
# 1015615 snRNA
#  362614 rRNA
#   10005 tRNA
#    3904 Simple_repeat
#    2356 SINE
#    2231 DNA
#     898 LTR
#     690 LINE
#     410 RNA
#     333 srpRNA
#     245 Low_complexity
#     101 Retroposon
#      51 Satellite
#      48 scRNA
#       5 LTR?
#       1 RC
#       1 DNA?
# SNRNP70_CTL_2.multimap.inrepeats.tsv
# 1159714 snRNA
#  744060 rRNA
#   13342 tRNA
#    6854 Simple_repeat
#    3047 SINE
#    2537 DNA
#    1157 RNA
#     952 LINE
#     951 srpRNA
#     849 LTR
#     310 Low_complexity
#     179 Retroposon
#      89 Satellite
#      51 scRNA
#       2 DNA?
#       1 LTR?

# how many U1 snRNA-mapping multimapping reads?
for tsv in *.multimap.inrepeats.tsv; do echo $tsv; grep snRNA $tsv | cut -f16 | grep -w "U1" | wc -l; done
# mVenus_CTL_1.multimap.inrepeats.tsv
# 1352
# mVenus_CTL_2.multimap.inrepeats.tsv
# 2568
# SNRNP70_CTL_1.multimap.inrepeats.tsv
# 1014969
# SNRNP70_CTL_2.multimap.inrepeats.tsv
# 1158252

# SNRNP70 CTL libraries mostly have multimapping reads, with more than half mapping to snRNAs (rep1: 68% rep2: 56%), followed by rRNAs (rep1: 24% rep2: 36%)
# importantly, SNRNP70 multimapping snRNA-mapped reads were predominated with U1 snRNA (rep1: 99% rep2: 99%)
# in other words, half of them in SNRNP70 are virtually U1 snRNA, while half are predominantly rRNAs among other repeat sequences.
# moreover, repeating from above, these multimapping reads are mostly found in repeats (rep1: 95%, rep2: 95%)
# taken together, SNRNP70 mostly binds to repeat regions during CTL, more than half of which are U1 snRNAs, followed by rRNAs.

# in contrast, mVenus multimappers were mostly rRNAs (rep1: 72%, rep2: 82%).
# these were followed either by a small proportion of snRNAs (rep1: 2%, rep2: 1%).

# note again, that there is a possibility of reads mapping to multiple regions if part of them maps to one and the other
# however it seems to be just for a very small number of reads, relatively:
for tsv in *.multimap.inrepeats.tsv; do echo $tsv; cut -f4 $tsv | sort | uniq -d | wc -l; done
# mVenus_CTL_1.multimap.inrepeats.tsv
# 1 (out of 53804)
# mVenus_CTL_2.multimap.inrepeats.tsv
# 1 (out of 153626)
# SNRNP70_CTL_1.multimap.inrepeats.tsv
# 2 (out of 1399508)
# SNRNP70_CTL_2.multimap.inrepeats.tsv
# 5 (out of 1934095)

#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=02_Cutadapt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%x_%A_%a.log
#SBATCH --mem=0
#SBATCH --time=100:00:00
#SBATCH --array=0-3

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate cutadapt
conda list

mkdir -p output

reads=(
SNRNP70-SA-rep1_L1_R1.fastq.gz 
SNRNP70-SA-rep2_L2_R1.fastq.gz
mVenus-SA-rep5_L1_R1.fastq.gz 
mVenus-SA-rep6_L2_R1.fastq.gz
)

FASTQ=input/${reads[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename ${FASTQ} .fastq.gz)

# 1st pass (-O 1)
cutadapt \
    -O 1 --match-read-wildcards --times 1 -e 0.1 --quality-cutoff 6 -m 18 \
    -o output/${SAMPLE}.fqTr.fq.gz \
    -a AGATCGGAAGAGCAC \
    -a GATCGGAAGAGCACA \
    -a ATCGGAAGAGCACAC \
    -a TCGGAAGAGCACACG \
    -a CGGAAGAGCACACGT \
    -a GGAAGAGCACACGTC \
    -a GAAGAGCACACGTCT \
    -a AAGAGCACACGTCTG \
    -a AGAGCACACGTCTGA \
    -a GAGCACACGTCTGAA \
    -a AGCACACGTCTGAAC \
    -a GCACACGTCTGAACT \
    -a CACACGTCTGAACTC \
    -a ACACGTCTGAACTCC \
    -a CACGTCTGAACTCCA \
    -a ACGTCTGAACTCCAG \
    -a CGTCTGAACTCCAGT \
    -a GTCTGAACTCCAGTC \
    -a TCTGAACTCCAGTCA \
    -a CTGAACTCCAGTCAC \
    ${FASTQ} \
    > output/${SAMPLE}.fqTr.metrics

# 2nd pass (-O 5)
cutadapt \
    -O 5 --match-read-wildcards --times 1 -e 0.1 --quality-cutoff 6 -m 18 \
    -o output/${SAMPLE}.trimmed.fq.gz \
    -a AGATCGGAAGAGCAC \
    -a GATCGGAAGAGCACA \
    -a ATCGGAAGAGCACAC \
    -a TCGGAAGAGCACACG \
    -a CGGAAGAGCACACGT \
    -a GGAAGAGCACACGTC \
    -a GAAGAGCACACGTCT \
    -a AAGAGCACACGTCTG \
    -a AGAGCACACGTCTGA \
    -a GAGCACACGTCTGAA \
    -a AGCACACGTCTGAAC \
    -a GCACACGTCTGAACT \
    -a CACACGTCTGAACTC \
    -a ACACGTCTGAACTCC \
    -a CACGTCTGAACTCCA \
    -a ACGTCTGAACTCCAG \
    -a CGTCTGAACTCCAGT \
    -a GTCTGAACTCCAGTC \
    -a TCTGAACTCCAGTCA \
    -a CTGAACTCCAGTCAC \
    output/${SAMPLE}.fqTr.fq.gz \
    > output/${SAMPLE}.fqTrTr.metrics

# CLIP-seq analysis
### Sample: SNRNP70
### Pipeline: Skipper (Boyle et al., Cell Genom., 2023 & Xu et al., STAR Protoc., 2024)
### Library prep: Denaturing Fluorescent Protein CLIP-seq
### Platform: DNA Nanoball Sequencing (BGI)
### Read type: Paired-end 100bp
### Date: 20251024

# Make mimic replicates
We only have one replicate per CLIP and input sample per treatment (e.g. non-stressed control, osmotic).
Since Skipper requires at least two replicates to run, but nonetheless generates some (not all) replicate-independent analysis, here I create mimic replicates, using 60% random reads from actual replicates.
Mimic replicates will be later labeled as "_2", while actual replicates for our analysis are "_1".
In the future, if for some reason the pipeline needs to be repeated, be careful **not to recreate** mimics using this script, so as to preserve reproducibility.

### Script: `01_MakeMimic.sh`
```bash
#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=Test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
conda activate skipper
ulimit -n 65536

# get random 60% reads to make mimic replicate
seqkit sample -p 0.6 input/N288_L1_1.fq.gz -o input/mimic_N288_L1_1.fq.gz
seqkit sample -p 0.6 input/N289_L1_1.fq.gz -o input/mimic_N289_L1_1.fq.gz
seqkit sample -p 0.6 input/N290_L1_1.fq.gz -o input/mimic_N290_L1_1.fq.gz
seqkit sample -p 0.6 input/N291_L1_1.fq.gz -o input/mimic_N291_L1_1.fq.gz
```

# Specify configuration
The manifest was configured as follows.
Note that for our experiment, our "size-matched inputs" are actually CLIP for mVenus HEK293T, while our samples are CLIP for SNRNP70-mVenus HEK293T.

### Manifest: `manifest.csv`
| Experiment     | Sample         | Cells   | Input_replicate | Input_adapter          | Input_fastq                  | CLIP_replicate | CLIP_adapter          | CLIP_fastq                   | Notes                          |
|----------------|----------------|---------|----------------|-----------------------|-------------------------------|----------------|----------------------|-------------------------------|--------------------------------|
| SNRNP70_CTL    | SNRNP70_CTL    | HEK293T | 1              | input/InvRiL19.fasta | input/N288_L1_1.fq.gz        | 1              | input/InvRiL19.fasta | input/N289_L1_1.fq.gz        |                                |
| SNRNP70_CTL    | SNRNP70_CTL    | HEK293T | 2              | input/InvRiL19.fasta | input/mimic_N288_L1_1.fq.gz  | 2              | input/InvRiL19.fasta | input/mimic_N289_L1_1.fq.gz  | mimic_60%_of_replicate_1       |
| SNRNP70_OSM    | SNRNP70_OSM    | HEK293T | 1              | input/InvRiL19.fasta | input/N290_L1_1.fq.gz        | 1              | input/InvRiL19.fasta | input/N291_L1_1.fq.gz        |                                |
| SNRNP70_OSM    | SNRNP70_OSM    | HEK293T | 2              | input/InvRiL19.fasta | input/mimic_N290_L1_1.fq.gz  | 2              | input/InvRiL19.fasta | input/mimic_N291_L1_1.fq.gz  | mimic_60%_of_replicate_1       |

We utilize the following fasta for adapaters in our experiment.

### Adapters: `InvRiL19.fasta`
```bash
>InvRiL19
AGATCGGAAGAGCACACGTC

```
Skipper's config file was customized as follows.
Note that for feature annotations, I utilized annotations for K562 cells, even though our data is from HEK293T.
Although both cell lines are human, note that this may lead to some differences.

### Config: `Skipper_config.py`
```bash
# Configuration file

import os
import sys
import glob

# Adjust REPO_PATH if necessary to make sure it contains the path to the pulled skipper repository when starting a skipper run
REPO_PATH = "."

########################################
# Customizable input files

MANIFEST = REPO_PATH + "/input/manifest.csv"
# Repeat table
REPEAT_TABLE = REPO_PATH + "/annotations/repeatmasker.grch38.tsv.gz"
# Genome fasta
GENOME = REPO_PATH + "/annotations/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
# STAR reference
STAR_DIR = REPO_PATH + "/annotations/genome_ref"
# Generated from STAR index
CHROM_SIZES = REPO_PATH + "/annotations/genome_ref/chrNameLength.txt"
# Use a GFF filtered for genes expressed in the cell type of interest
GFF = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.gff3.gz"
# Customizable, with defaults
BLACKLIST = REPO_PATH + "/annotations/encode3_eclip_blacklist.bed" # set to None for no blacklisting
GENE_SETS = REPO_PATH + "/annotations/c5.go.v7.5.1.symbols.gmt"
GENE_SET_REFERENCE = REPO_PATH + "/annotations/encode3_go_terms.reference.tsv.gz"
GENE_SET_DISTANCE = REPO_PATH + "/annotations/encode3_go_terms.jaccard_index.rds"
# Ranked list of gene and transcript types found in GFF annotations
ACCESSION_RANKINGS = REPO_PATH + "/annotations/accession_type_ranking.txt"

########################################
# Customizable parameters

# Information about CLIP library
UMI_SIZE = 10
# Single-end: enter 1. Paired-end: enter read (1 or 2) corresponding to crosslink site
INFORMATIVE_READ = 1
# Internal use
UNINFORMATIVE_READ = 3 - INFORMATIVE_READ
# Use multiple input replicates to estimate overdispersion (preferred), or use multiple CLIP replicates
# Skipper requires replicates to model the variance in read counts.
OVERDISPERSION_MODE = "input" # input or clip

########################################
# Intermediate files and scripts for the Skipper run, user setup not required. Adjust the paths if necessary.

# Skipper will partition the transcriptome and create feature annotations from the GFF
PARTITION = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.bed.gz"
FEATURE_ANNOTATIONS = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.features.tsv.gz"
# Skipper will sort the repeat table.
REPEAT_BED = REPO_PATH + "/annotations/repeatmasker.grch38.sort.unique.bed.gz"
# The directory contains Skipper scripts
TOOL_DIR = REPO_PATH + "/tools"

```

# Check data integrity
Since the raw data has been copied, pasted, moved, etc., I checked the md5.
The values were the same as the ones when we received them from BGI.
The results of the check can be found at `md5CheckSums.log`.

### Script: `02_CheckData.sh`
```bash
#!/bin/bash

LOGFILE="./input/md5CheckSums.log"
echo "MD5 Checksum checked on: $(date)" > "$LOGFILE"
echo "----------------------------------------" >> "$LOGFILE"
md5sum ./input/*.fq.gz >> "$LOGFILE"
echo "MD5 checksums written to $LOGFILE"
```

# Test Skipper
As recommended by the documentation, test Skipper first for any conflicts (e.g. directories, unique naming of replicates, etc.)

### Script: `03_TestSkipper.sh`
```bash
#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=Test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%j.log
#SBATCH --mem=0
#SBATCH --time=100:00:00

source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
conda activate skipper
ulimit -n 65536

snakemake -ns Skipper.py -j 1
```
# Switch to interactive slurm cluster
Although Skipper's documentation recommended adding a slurm cluster specification command in addition to the Skipper command, I found that it somehow kept creating errors for the `OSM` samples, especially when executing the rule `partition_bam_reads` and `quantify_repeats`; when I execute the shell commands from those rules manually, however, the outputs were created fine and the pipeline could proceed. Therefore, better to instead execute it as below. From our computational server's development node, I switch to an interactive node, assigned to me by slurm.

### Script: `04_ActivateNode.sh`
```bash
#!/bin/bash

srun -p cluster_intr -c 52 --time=10:00:00 --pty bash -l
```

# Run Skipper
I ran the Skipper pipeline as follows.

### Script: `05_RunSkipper.sh`
```bash
#!/bin/bash

source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
conda activate skipper
snakemake -kps Skipper.py -j 30 -w 30
```

# Rename outputs
Skipper's documentation recommends renaming the resultant `/output` directory.
This will avoid overwriting the results, for future pipeline executions.
As personal preference, I also renamed some other files as specified below.
The files can be viewed at `/project/okamura-lab-imaging/Patrick_forPatrick/Patrick_Bioinformatics/SNRNP70`.
Moreover, the log `2025-10-24T132715.983234.snakemake.log` from the hidden `.snakemake` can also be viewed there.
Additionally, I manually copy-pasted my actual log into `20251024_ManualLog.log` (I later realized `stderr` and `stdout` becomes empty). 

```bash
mv input 20251024_input
mv output 20251024_output
mv benchmarks 20251024_benchmarks
mv Skipper.py 20251024_Skipper.py
mv Skipper_config.py 20251024_Skipper_config.py
```

### Patrick, (last edited 20251103)

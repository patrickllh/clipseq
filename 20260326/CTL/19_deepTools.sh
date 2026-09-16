# #!/bin/bash

# #SBATCH --partition=cluster_long
# #SBATCH --job-name=19_deepTools
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=52
# #SBATCH --output=%x_%j.log
# #SBATCH --mem=0
# #SBATCH --time=100:00:00

# set -euo pipefail
# source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
# ulimit -n 65536
# conda activate deeptools
# conda list

# mkdir -p output

# awk '$0 ~ /^#/ \
# || ($3=="transcript" && $0 ~ /transcript_type "protein_coding"/) \
# || ($3=="exon" && $0 ~ /transcript_type "protein_coding"/)' \
# input/gencode.v38.annotation.gtf > output/mrna.gtf

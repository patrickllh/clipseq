#!/bin/bash
set -euo pipefail
mkdir -p inputFolder
htseq-clip count -i SNRNP70_OSM_IP_1.crosslinks.gz -a slidingwindows.gz -o inputFolder/SNRNP70_OSM_IP_1_counts.gz
htseq-clip count -i SNRNP70_OSM_IP_2.crosslinks.gz -a slidingwindows.gz -o inputFolder/SNRNP70_OSM_IP_2_counts.gz
htseq-clip count -i SNRNP70_OSM_IN_1.crosslinks.gz -a slidingwindows.gz -o inputFolder/SNRNP70_OSM_IN_1_counts.gz
htseq-clip count -i SNRNP70_OSM_IN_2.crosslinks.gz -a slidingwindows.gz -o inputFolder/SNRNP70_OSM_IN_2_counts.gz

#!/bin/bash
set -euo pipefail

htseq-clip extract -e 1 -s s --ignore -c 52 -i SNRNP70_OSM_IP_1.genome.Aligned.sort.dedup.bam -o SNRNP70_OSM_IP_1.crosslinks.gz
htseq-clip extract -e 1 -s s --ignore -c 52 -i SNRNP70_OSM_IP_2.genome.Aligned.sort.dedup.bam -o SNRNP70_OSM_IP_2.crosslinks.gz
htseq-clip extract -e 1 -s s --ignore -c 52 -i SNRNP70_OSM_IN_1.genome.Aligned.sort.dedup.bam -o SNRNP70_OSM_IN_1.crosslinks.gz
htseq-clip extract -e 1 -s s --ignore -c 52 -i SNRNP70_OSM_IN_2.genome.Aligned.sort.dedup.bam -o SNRNP70_OSM_IN_2.crosslinks.gz

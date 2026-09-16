#!/bin/bash
set -euo pipefail

htseq-clip createMatrix -i inputFolder -b SNRNP70 -o matrix.gz

#!/bin/bash

htseq-clip annotation -g gencode.v38.annotation.gff3.gz --splitExons --unsorted -o annotation.gz

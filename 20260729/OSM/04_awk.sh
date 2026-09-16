#!/bin/bash

for f in output/*peaks.normed.compressed.bed; do
    awk -F'\t' '$5 > 0' "$f" > "${f%.bed}.log2FCpos.bed"
done

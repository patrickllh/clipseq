#!/bin/bash

LOGFILE="./input/md5CheckSums.log"
echo "MD5 Checksum checked on: $(date)" > "$LOGFILE"
echo "----------------------------------------" >> "$LOGFILE"
md5sum ./input/*.fq.gz >> "$LOGFILE"
echo "MD5 checksums written to $LOGFILE"

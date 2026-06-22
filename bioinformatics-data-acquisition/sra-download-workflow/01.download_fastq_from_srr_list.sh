#!/bin/bash

# Usage:
# bash 01.download_fastq_from_srr_list.sh SRR_Acc_List.txt

if [ -z "$1" ]; then
    echo "Usage: bash 01.download_fastq_from_srr_list.sh SRR_Acc_List.txt"
    exit 1
fi

srr_list=$1

mkdir -p fastq

while read id
do
    [ -z "$id" ] && continue

    echo "Downloading ${id}..."

    prefetch "$id"

    fasterq-dump "$id" \
        --split-files \
        --threads 4 \
        -O fastq

    echo "Finished ${id}"

done < "$srr_list"

echo "Done."

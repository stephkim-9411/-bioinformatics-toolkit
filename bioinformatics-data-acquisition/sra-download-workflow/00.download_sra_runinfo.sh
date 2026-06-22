#!/bin/bash

# Usage:
# bash 00.download_sra_runinfo.sh PRJNA1051047

if [ -z "$1" ]; then
    echo "Usage: bash 00.download_sra_runinfo.sh <BioProject_ID>"
    exit 1
fi

project_id=$1

echo "Downloading RunInfo for ${project_id}..."

curl -L \
  -o runinfo.csv \
  "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${project_id}"

if [ ! -s runinfo.csv ]; then
    echo "Error: runinfo.csv was not downloaded."
    exit 1
fi

echo "Extracting SRR list..."

cut -d',' -f1 runinfo.csv | grep SRR > SRR_Acc_List.txt

echo "Done."
echo "Generated:"
echo "  runinfo.csv"
echo "  SRR_Acc_List.txt"

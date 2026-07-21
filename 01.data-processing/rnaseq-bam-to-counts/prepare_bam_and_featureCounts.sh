#!/bin/bash

set -euo pipefail

# ============================================================
# BAM quality check, sorting, indexing, and gene quantification
#
# Usage:
# bash prepare_bam_and_featureCounts.sh \
#     /path/to/input_bams \
#     /path/to/annotation.gtf \
#     /path/to/output \
#     16
# ============================================================

INPUT_DIR=$1
GTF=$2
OUTPUT_DIR=$3
THREADS=${4:-16}

SORTED_DIR="${OUTPUT_DIR}/sorted_bam"
COUNTS_DIR="${OUTPUT_DIR}/featureCounts_results"

mkdir -p "$SORTED_DIR" "$COUNTS_DIR"

#----------------------------------------------------
# Check tools
#----------------------------------------------------

command -v samtools >/dev/null || { echo "samtools not found"; exit 1; }
command -v featureCounts >/dev/null || { echo "featureCounts not found"; exit 1; }

#----------------------------------------------------
# Check chromosome naming
#----------------------------------------------------

GTF_CHR=$(awk '!/^#/ {print $1; exit}' "$GTF")

if [[ "$GTF_CHR" == chr* ]]; then
    GTF_STYLE="chr"
else
    GTF_STYLE="no_chr"
fi

FIRST_BAM=$(ls "$INPUT_DIR"/*.bam | head -1)

BAM_CHR=$(samtools view -H "$FIRST_BAM" |
    awk '/^@SQ/{
        sub("SN:","",$2)
        print $2
        exit
    }')

if [[ "$BAM_CHR" == chr* ]]; then
    BAM_STYLE="chr"
else
    BAM_STYLE="no_chr"
fi

USE_ALIAS=false

if [ "$BAM_STYLE" != "$GTF_STYLE" ]; then

    echo "Chromosome naming differs."
    echo "Creating alias file..."

    ALIAS_FILE="${OUTPUT_DIR}/chromosome_alias.txt"

    for i in {1..22}; do
        echo "chr${i},${i}"
    done > "$ALIAS_FILE"

    echo "chrX,X" >> "$ALIAS_FILE"
    echo "chrY,Y" >> "$ALIAS_FILE"
    echo "chrM,MT" >> "$ALIAS_FILE"

    USE_ALIAS=true

fi

echo "Sorting and indexing BAM files..."

BAMS=()

for BAM in "$INPUT_DIR"/*.bam
do

    SAMPLE=$(basename "$BAM" .bam)
    SORTED_BAM="${SORTED_DIR}/${SAMPLE}.sorted.bam"

    samtools quickcheck "$BAM"

    samtools sort \
        -@ "$THREADS" \
        -o "$SORTED_BAM" \
        "$BAM"

    samtools index "$SORTED_BAM"

    BAMS+=("$SORTED_BAM")

done

echo "Running featureCounts..."

FC_OPTIONS=(
    -T "$THREADS"
    -p
    --countReadPairs
    -t exon
    -g gene_id
    -a "$GTF"
    -o "${COUNTS_DIR}/gene_counts.txt"
)

if [ "$USE_ALIAS" = true ]; then
    FC_OPTIONS+=(-A "$ALIAS_FILE")
fi

featureCounts \
    "${FC_OPTIONS[@]}" \
    "${BAMS[@]}"

echo
echo "Finished!"
echo

grep "^Assigned" "${COUNTS_DIR}/gene_counts.txt.summary"

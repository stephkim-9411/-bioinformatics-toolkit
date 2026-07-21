## Post-alignment Processing (BAM → Gene Counts)


This workflow starts from aligned RNA-seq BAM files.

### Steps

1. Verify BAM files exist.
2. Check BAM integrity using `samtools quickcheck`.
3. Sort BAM files by genomic coordinates.
4. Index sorted BAM files.
5. Quantify gene-level read counts using `featureCounts`.

### Required Tools

- **samtools** (v1.10 or later recommended)
  - BAM integrity check
  - BAM sorting
  - BAM indexing

- **Subread** (includes `featureCounts`)
  - Gene-level read quantification

### Input

- Aligned RNA-seq BAM files (`*.bam`)
- GENCODE GTF annotation (`*.gtf`)

### Output

```text
sorted_bam/
├── sample1.sorted.bam
├── sample1.sorted.bam.bai
└── ...

featureCounts_results/
├── gene_counts.txt
└── gene_counts.txt.summary

logs/
└── featureCounts.log
```

### Notes

- Ensure the BAM files and GTF annotation use the **same reference genome** (e.g., GRCh38).
- Ensure chromosome naming is consistent between BAM and GTF (e.g., `chr1` vs `1`).
- This workflow is designed for **paired-end RNA-seq** data.

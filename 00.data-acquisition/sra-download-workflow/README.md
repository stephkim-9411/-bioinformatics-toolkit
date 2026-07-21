# SRA Data Download Workflow

Simple scripts for downloading public sequencing datasets from NCBI SRA.

---

## Overview

This workflow allows you to:

1. Download metadata from an NCBI BioProject
2. Extract all SRR accession IDs
3. Download FASTQ files for all samples

Workflow:

```text
BioProject ID
      │
      ▼
runinfo.csv
      │
      ▼
SRR_Acc_List.txt
      │
      ▼
FASTQ files
```

---

## Requirements

### Install Miniconda

If Conda is not installed:

https://docs.conda.io/en/latest/miniconda.html

---

### Create Environment

Install NCBI SRA Toolkit:

```bash
conda create -n sra-tools -c bioconda sra-tools -y
```

Activate environment:

```bash
conda activate sra-tools
```

---

## Step 1. Download RunInfo and Generate SRR List

Run:

```bash
bash 00.download_sra_runinfo.sh {project_id}
```

Generated files:

```text
runinfo.csv
SRR_Acc_List.txt
```

### What are these files?

#### runinfo.csv

Contains metadata for every sequencing run.

#### SRR_Acc_List.txt

Contains only SRR accession IDs.

This file will be used in the next step.

---

## Step 2. Download FASTQ Files

Use the SRR list generated in Step 1.

Run:

```bash
bash 01.download_fastq_from_srr_list.sh SRR_Acc_List.txt
```

Important:

❗ Do NOT use `runinfo.csv`

✅ Use:

```bash
SRR_Acc_List.txt
```

The script reads each SRR accession and downloads the corresponding FASTQ files.

---

## Output

After completion, FASTQ files will be saved in the `fastq/` directory.

## Author

**Stephanie Kim**  
Ph.D. Student, Human Genetics and Genomics  
University of Miami

## Citation

If you use this pipeline in your research, please cite this repository:

> Kim, S. (2026). *SRA Data Download Pipeline*. GitHub.  
> https://github.com/stephkim-9411/bioinformatics-toolkit

## Support

If you find this repository helpful, please consider giving it a ⭐.


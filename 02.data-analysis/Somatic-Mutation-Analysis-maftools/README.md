# maftools tutorial: MAF preparation, mutation analysis, and optional CNV

A step-by-step R Markdown tutorial using the **public TCGA AML example supplied by maftools**. No private input files are needed for the default workflow. The example is real public research data, not simulated data.


This tutorial uses the **public TCGA AML example supplied by maftools** to explore somatic mutation patterns across samples. You will learn how to identify frequently mutated genes, visualize mutation types and protein positions, and summarize alterations across a cohort.

Optional sections demonstrate how to include gene-level copy-number alterations and explore previously generated OncoKB annotations using your own data.

## What can you explore?

| Analysis | What you can see |
|---|---|
| **Oncoplot** | Which genes are frequently altered, which samples carry those alterations, and which alteration types occur |
| **Selected-gene oncoplot** | Mutation patterns in a gene panel of your choice |
| **Mutation summary** | The distribution of mutation consequences and variant types, along with mutation counts per sample |
| **Gene and sample summaries** | How many samples have alterations in each gene and how many mutations are recorded in each sample |
| **Protein lollipop plots** | Where mutations fall along a protein, including recurrent positions and their relationship to annotated domains |
| **Rainfall plots** | How mutations are spaced along the genome and where closely spaced mutations occur |
| **Optional CNV integration** | Amplifications and deletions alongside sequence mutations in the same oncoplot |
| **Optional OncoKB summaries** | Oncogenicity, functional effects, and reported therapeutic evidence for previously annotated variants |
| **Optional mutational signatures** | Patterns of nucleotide substitutions and their similarity to reference signatures |

The example dataset contains real public research data, not simulated data. CNV and OncoKB analyses require additional inputs and are disabled by default.

Mutation counts are not equivalent to tumor mutational burden (TMB). Rainfall clusters and signature similarities are exploratory findings and do not establish a biological mechanism on their own.

## Files

| File | Purpose |
|---|---|
| `01_MAF_preparation_and_annotation.md` | VCF-to-MAF preparation, merging per-sample MAFs, and optional OncoKB annotation |
| `02_maftools_analysis.R` | R script for mutation analysis, visualization, and optional CNV integration |
| `data/tcga_laml.maf.gz` | Public maftools example MAF, if included; otherwise loaded from the installed package |
| `data/README.md` | Documentation for bundled example data |
| `data/maftools_LICENSE.txt` | Upstream license notice for bundled example data |

## Run the example

1. Download this folder and open `02_maftools_analysis.R` in RStudio.
2. Set your working directory to the folder containing the script.
3. Install the required packages once:

   ```r
   install.packages(c("BiocManager", "data.table", "dplyr", "ggplot2", "R.utils"))
   BiocManager::install("maftools")

The tutorial reads the included example file. If it is missing, it falls back to the copy installed with maftools. A different installed package version may contain a different example file.

## Analyses

- VCF-to-MAF preparation instructions and combining per-sample MAFs.
- Top-gene and selected-gene oncoplots.
- Mutation summary dashboard and gene/sample summary tables.
- Top-gene protein lollipop plots.
- Sample-level rainfall plots.
- Optional CNV Amp/Del integration.
- Optional OncoKB evidence/effect plots and variant tables.
- Optional SBS signature analysis.

## Expected results

Default outputs include `Oncoplot.pdf`, `Selected_genes_oncoplot.pdf`, `Mutation_MAF_summary.pdf`, gene/sample CSVs, `Top20_LollipopPlots.pdf`, `RainfallPlots.pdf`, `Plot_error_log.csv`, `MAF_objects.rds`, and `sessionInfo.txt`.

The MAF contains 2,207 input rows from 193 sample IDs. maftools filters consequences and deduplicates records, so output mutation counts can differ. Lollipop and rainfall failures are logged. No pre-rendered results are included: R was unavailable in the preparation environment, so this package has not been executed or knitted there.

## Use your own data

Set `use_example <- FALSE`, update `maf_file` and `output_dir`, and set the correct genome build. Enable `use_cnv` only with matched gene-level CNV calls, and `use_oncokb` only with an OncoKB-annotated MAF. The AML example has neither CNV calls nor OncoKB annotation, so those sections remain disabled in the example run.

## Publish on GitHub

Upload the contents of this folder to your chosen tutorial directory, preserving `data/`. The included public example should remain separate from private analysis inputs. Run and inspect the tutorial locally before adding a rendered HTML or figures as expected results.

## Acknowledgments

Data and software: [maftools](https://github.com/PoisonAlien/maftools). Use `citation("maftools")` for the software citation. Data provenance and the upstream license notice are in `data/`.

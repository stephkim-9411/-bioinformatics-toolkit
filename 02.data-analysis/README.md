# RNA-seq Differential Expression Tutorial

A beginner-friendly, reproducible R workflow for **Control vs Case RNA-seq differential expression** using both **edgeR** and **DESeq2**, followed by gene annotation, GO, KEGG, STRING, and GSEA.

The repository intentionally separates differential expression from organism-specific annotation. This allows the DEG scripts to operate on raw gene identifiers first and lets the user select the organism later.

## Workflow

```text
Raw count matrix
      |
      v
01. Prepare metadata
      |
      +-------------------+
      |                   |
      v                   v
02A. edgeR           02B. DESeq2
      |                   |
      +---------+---------+
                |
                v
       02C. Compare methods
                |
                v
       03. Gene annotation
       human / mouse / rat
                |
        +-------+-------+
        |       |       |
        v       v       v
      04 GO   05 KEGG  06 STRING
                |
                v
             07 GSEA
```

## Repository structure

```text
RNAseq_DEG_Tutorial/
├── README.md
├── data/
│   ├── example_counts.csv
│   └── example_metadata.csv
├── scripts/
│   ├── 00_install_packages.R
│   ├── 01_prepare_metadata.R
│   ├── 02A_edgeR_DEG.R
│   ├── 02B_DESeq2_DEG.R
│   ├── 02C_compare_DEG_methods.R
│   ├── functions/
│   │   └── volcano_plot.R
│   ├── 03_gene_annotation.R
│   ├── 04_GO_enrichment.R
│   ├── 05_KEGG_enrichment.R
│   ├── 06_STRING_network.R
│   └── 07_GSEA.R
└── results/
    ├── DEG/
    ├── annotation/
    ├── GO/
    ├── KEGG/
    ├── STRING/
    └── GSEA/
```

## Example data

`example_counts.csv` contains **synthetic RNA-seq counts** for 100 real mouse gene symbols. The expression values are simulated and should not be interpreted biologically.

Samples are named:

```text
Control_1 ... Control_6
Case_1    ... Case_6
```

This naming convention allows `01_prepare_metadata.R` to automatically create the group information.

## Before starting

Open R/RStudio with the repository root as the working directory.

Install packages once:

```r
source("scripts/00_install_packages.R")
```

## 1. Prepare metadata

```r
source("scripts/01_prepare_metadata.R")
```

The script extracts sample columns from the count matrix and identifies `Control` and `Case` from the sample names.

For your own data, edit the pattern-matching rules if your samples use a different naming convention.

## 2A. edgeR

```r
source("scripts/02A_edgeR_DEG.R")
```

This tutorial uses the edgeR quasi-likelihood workflow:

1. `DGEList`
2. `filterByExpr`
3. TMM normalization
4. dispersion estimation
5. `glmQLFit`
6. `glmQLFTest`

Positive `logFC` means higher expression in **Case** relative to **Control**.

## 2B. DESeq2

```r
source("scripts/02B_DESeq2_DEG.R")
```

Positive `log2FoldChange` means higher expression in **Case** relative to **Control**.

Both edgeR and DESeq2 require **raw integer counts**, not TPM, FPKM, CPM, or previously normalized expression values.

## Change the fold-change threshold

At the top of both DEG scripts:

```r
pvalue_threshold <- 0.05
logfc_threshold <- 1
```

For a stricter fold-change cutoff:

```r
logfc_threshold <- 2
```

The same value controls the DEG classification and volcano-plot vertical threshold.

## Volcano plots

The shared function is stored in:

```text
scripts/functions/volcano_plot.R
```

The plot follows this convention:

- red: up-regulated
- blue: down-regulated
- black: non-significant
- y-axis: `-log2(P-value)`
- dashed horizontal line: selected P-value threshold
- dashed vertical lines: selected absolute log2 fold-change threshold
- labels: top positive/negative genes plus genes passing the selected thresholds

To reproduce a fixed y-axis similar to the original analysis:

```r
volcano_y_max <- 15
```

To let R determine the range automatically:

```r
volcano_y_max <- NULL
```

## 2C. Compare edgeR and DESeq2

Run both DEG scripts first, then:

```r
source("scripts/02C_compare_DEG_methods.R")
```

The script compares fold-change estimates and creates a table of genes significant in both methods.

Agreement between methods is informative, but a gene does not need to be significant in both methods to be biologically valid. edgeR and DESeq2 use different statistical models and estimation procedures.

## 3. Select organism and annotate genes

Annotation is deliberately performed **after DEG analysis**.

At the top of `03_gene_annotation.R`:

```r
organism <- "mouse"
gene_id_type <- "SYMBOL"
```

Supported tutorial organisms:

```r
organism <- "human"
organism <- "mouse"
organism <- "rat"
```

Common identifier types include:

```r
gene_id_type <- "SYMBOL"
gene_id_type <- "ENSEMBL"
gene_id_type <- "ENTREZID"
```

The organism determines the annotation resource and downstream identifiers:

| Organism | Annotation database | KEGG | STRING taxonomy |
|---|---|---|---:|
| Human | `org.Hs.eg.db` | `hsa` | 9606 |
| Mouse | `org.Mm.eg.db` | `mmu` | 10090 |
| Rat | `org.Rn.eg.db` | `rno` | 10116 |

The original gene identifier is retained rather than replaced.

Run:

```r
source("scripts/03_gene_annotation.R")
```

## 4. GO enrichment

```r
source("scripts/04_GO_enrichment.R")
```

This script performs **over-representation analysis (ORA)** using significant DEGs and GO Biological Process terms.

## 5. KEGG enrichment

```r
source("scripts/05_KEGG_enrichment.R")
```

This is also an **over-representation analysis**. It should not be confused with GSEA.

## 6. STRING network

```r
source("scripts/06_STRING_network.R")
```

The script maps significant gene symbols to STRING and creates a protein-protein interaction network.

STRING requires internet access when retrieving network information.

## 7. GSEA

```r
source("scripts/07_GSEA.R")
```

Unlike GO/KEGG ORA, GSEA uses the **full ranked gene list**, not only genes passing a significance threshold.

The tutorial ranks genes by DESeq2 `log2FoldChange`.

## Recommended run order

```r
source("scripts/00_install_packages.R")
source("scripts/01_prepare_metadata.R")

source("scripts/02A_edgeR_DEG.R")
source("scripts/02B_DESeq2_DEG.R")
source("scripts/02C_compare_DEG_methods.R")

source("scripts/03_gene_annotation.R")
source("scripts/04_GO_enrichment.R")
source("scripts/05_KEGG_enrichment.R")
source("scripts/06_STRING_network.R")
source("scripts/07_GSEA.R")
```

## Important notes

### Raw counts

edgeR and DESeq2 should receive raw integer read counts.

Do **not** use:

```text
TPM
FPKM
RPKM
log-transformed counts
CPM
```

as the count matrix for these DEG models.

### P-value versus adjusted P-value

The tutorial uses adjusted P-values/FDR to define DEGs. The volcano visualization retains the original analysis style and displays raw P-values on the y-axis.

### Synthetic example data

The included counts are simulated for teaching. Although real mouse gene symbols are used so annotation tools can run, enrichment results from this example should **not** be interpreted as biological findings.

### Organism consistency

After changing the organism in the annotation step, use the same organism in GO, KEGG, STRING, and GSEA scripts.

## Extending the tutorial

Natural next modules include:

- PCA and sample QC
- heatmaps
- batch-effect covariates
- paired designs
- sex-by-condition interaction models
- multiple experimental groups
- pathway visualization
- enrichment using up- and down-regulated genes separately

## License / educational use

This repository is intended as an educational template. Users should adapt the experimental design, filtering strategy, statistical thresholds, organism, and annotation resources to their own study.

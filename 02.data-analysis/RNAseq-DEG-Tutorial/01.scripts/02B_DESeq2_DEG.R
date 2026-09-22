# ============================================================
# 02B. Differential Expression Analysis Using DESeq2
# Identifies genes differentially expressed between two groups.
# ============================================================

library(data.table)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Load the custom volcano plot function
source("scripts/functions/volcano_plot.R")

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

counts_file <- "data/example_counts.csv"
metadata_file <- "data/example_metadata.csv"

gene_column <- "Gene_ID"
reference_group <- "Control"
comparison_group <- "Case"

pvalue_threshold <- 0.05
logfc_threshold <- 1  # Change to 2 for |log2FC| > 2
volcano_y_max <- 70   # Set to NULL for an automatic y-axis

# ------------------------------------------------------------
# 2. Read count and metadata files
# ------------------------------------------------------------

counts_raw <- data.table::fread(counts_file)
metadata <- data.table::fread(metadata_file)

# ------------------------------------------------------------
# 3. Prepare the count matrix
# ------------------------------------------------------------

# Save the gene IDs
gene_ids <- counts_raw[[gene_column]]

# Remove the gene ID column and convert the counts to a matrix
counts <- as.matrix(
  counts_raw[
    ,
    setdiff(names(counts_raw), gene_column),
    with = FALSE
  ]
)

# Assign gene IDs as row names
rownames(counts) <- gene_ids

# DESeq2 requires integer counts
storage.mode(counts) <- "integer"

# ------------------------------------------------------------
# 4. Check sample names and order
# ------------------------------------------------------------

if (!identical(colnames(counts), metadata$Sample_ID)) {
  stop(
    "Sample order/names in counts do not exactly match metadata$Sample_ID."
  )
}

# ------------------------------------------------------------
# 5. Prepare the metadata
# ------------------------------------------------------------

metadata <- as.data.frame(metadata)

# Assign sample IDs as metadata row names
rownames(metadata) <- metadata$Sample_ID

# Set Control as the reference and Case as the comparison
metadata$Group <- factor(
  metadata$Group,
  levels = c(reference_group, comparison_group)
)

# ------------------------------------------------------------
# 6. Create the DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Group
)

# ------------------------------------------------------------
# 7. Filter genes with very low counts
# Keep genes with at least 10 total reads across all samples.
# ------------------------------------------------------------

keep <- rowSums(DESeq2::counts(dds)) >= 10

dds <- dds[keep, ]

message(
  "Number of genes retained after filtering: ",
  nrow(dds)
)

# ------------------------------------------------------------
# 8. Run differential expression analysis
# ------------------------------------------------------------

dds <- DESeq2::DESeq(dds)

# Case is compared against Control
res <- DESeq2::results(
  dds,
  contrast = c(
    "Group",
    comparison_group,
    reference_group
  ),
  alpha = pvalue_threshold
)

# ------------------------------------------------------------
# 9. Convert results to a data frame
# ------------------------------------------------------------

DESeq2_results <- as.data.frame(res) |>
  tibble::rownames_to_column(
    var = gene_column
  )

# ------------------------------------------------------------
# 10. Classify gene regulation
# ------------------------------------------------------------

DESeq2_results <- DESeq2_results |>
  dplyr::mutate(
    Regulation = dplyr::case_when(
      !is.na(padj) &
        padj < pvalue_threshold &
        log2FoldChange > logfc_threshold ~ "Up",
      
      !is.na(padj) &
        padj < pvalue_threshold &
        log2FoldChange < -logfc_threshold ~ "Down",
      
      TRUE ~ "Not significant"
    )
  )

# ------------------------------------------------------------
# 11. Create the results directory
# ------------------------------------------------------------

dir.create(
  "results/DEG",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 12. Save the complete DESeq2 results
# ------------------------------------------------------------

write.csv(
  DESeq2_results,
  file = "results/DEG/DESeq2_Case_vs_Control.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Create the volcano plot
# Uses adjusted p-values (padj).
# ------------------------------------------------------------

p_DESeq2 <- create_volcano_plot(
  DESeq2_results,
  gene_column = gene_column,
  logfc_column = "log2FoldChange",
  pvalue_column = "padj",
  logfc_threshold = logfc_threshold,
  pvalue_threshold = pvalue_threshold,
  y_max = volcano_y_max,
  title = "Volcano Plot: DESeq2 Case vs Control"
)

print(p_DESeq2)

# ------------------------------------------------------------
# 14. Save the volcano plot
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = "results/DEG/DESeq2_volcano.png",
  plot = p_DESeq2,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Count genes in each regulation category
# ------------------------------------------------------------

regulation_summary <- table(
  DESeq2_results$Regulation
)

print(regulation_summary)

# ------------------------------------------------------------
# 16. Save the regulation summary
# ------------------------------------------------------------

write.csv(
  as.data.frame(regulation_summary),
  file = "results/DEG/DESeq2_regulation_summary.csv",
  row.names = FALSE
)

message(
  "DESeq2 differential expression analysis completed."
)

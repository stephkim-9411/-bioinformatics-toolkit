# ============================================================
# 02A. Differential Expression Analysis Using edgeR
# Identifies genes differentially expressed between two groups.
# ============================================================

# Set the project directory

library(data.table)
library(dplyr)
library(edgeR)
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
volcano_y_max <- 40   # Set to NULL for an automatic y-axis

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

# Remove the gene ID column and convert counts to a matrix
counts <- as.matrix(
  counts_raw[
    ,
    setdiff(names(counts_raw), gene_column),
    with = FALSE
  ]
)

# Assign gene IDs as row names
rownames(counts) <- gene_ids

# edgeR requires integer count values
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
# 5. Prepare the sample groups
# ------------------------------------------------------------

metadata$Group <- factor(
  metadata$Group,
  levels = c(reference_group, comparison_group)
)

# ------------------------------------------------------------
# 6. Create the edgeR DGEList object
# ------------------------------------------------------------

dge <- edgeR::DGEList(
  counts = counts,
  group = metadata$Group
)

# ------------------------------------------------------------
# 7. Filter genes with low expression
# ------------------------------------------------------------

keep <- edgeR::filterByExpr(
  dge,
  group = metadata$Group
)

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

message(
  "Number of genes retained after filtering: ",
  nrow(dge)
)

# ------------------------------------------------------------
# 8. Calculate TMM normalization factors
# ------------------------------------------------------------

dge <- edgeR::normLibSizes(
  dge,
  method = "TMM"
)

# ------------------------------------------------------------
# 9. Create the design matrix
# ------------------------------------------------------------

design <- stats::model.matrix(
  ~ Group,
  data = metadata
)

print(design)

# ------------------------------------------------------------
# 10. Estimate dispersion
# ------------------------------------------------------------

dge <- edgeR::estimateDisp(
  dge,
  design
)

# ------------------------------------------------------------
# 11. Fit the quasi-likelihood model
# ------------------------------------------------------------

fit <- edgeR::glmQLFit(
  dge,
  design,
  robust = TRUE
)

# ------------------------------------------------------------
# 12. Test Case versus Control
# ------------------------------------------------------------

test <- edgeR::glmQLFTest(
  fit,
  coef = paste0("Group", comparison_group)
)

# ------------------------------------------------------------
# 13. Extract differential expression results
# ------------------------------------------------------------

edgeR_results <- edgeR::topTags(
  test,
  n = Inf
)$table |>
  tibble::rownames_to_column(
    var = gene_column
  )

# ------------------------------------------------------------
# 14. Classify gene regulation
# ------------------------------------------------------------

edgeR_results <- edgeR_results |>
  dplyr::mutate(
    Regulation = dplyr::case_when(
      !is.na(FDR) &
        FDR < pvalue_threshold &
        logFC > logfc_threshold ~ "Up",
      
      !is.na(FDR) &
        FDR < pvalue_threshold &
        logFC < -logfc_threshold ~ "Down",
      
      TRUE ~ "Not significant"
    )
  )

# ------------------------------------------------------------
# 15. Create the results directory
# ------------------------------------------------------------

dir.create(
  "results/DEG",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 16. Save the complete edgeR results
# ------------------------------------------------------------

write.csv(
  edgeR_results,
  file = "results/DEG/edgeR_Case_vs_Control.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 17. Create the volcano plot
# Uses adjusted p-values (FDR).
# ------------------------------------------------------------

p_edgeR <- create_volcano_plot(
  edgeR_results,
  gene_column = gene_column,
  logfc_column = "logFC",
  pvalue_column = "FDR",
  logfc_threshold = logfc_threshold,
  pvalue_threshold = pvalue_threshold,
  y_max = volcano_y_max,
  title = "Volcano Plot: edgeR Case vs Control"
)

print(p_edgeR)

# ------------------------------------------------------------
# 18. Save the volcano plot
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = "results/DEG/edgeR_volcano.png",
  plot = p_edgeR,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 19. Count genes in each regulation category
# ------------------------------------------------------------

regulation_summary <- table(
  edgeR_results$Regulation
)

print(regulation_summary)

# ------------------------------------------------------------
# 20. Save the regulation summary
# ------------------------------------------------------------

write.csv(
  as.data.frame(regulation_summary),
  file = "results/DEG/edgeR_regulation_summary.csv",
  row.names = FALSE
)

message(
  "edgeR differential expression analysis completed."
)
if (!identical(colnames(counts), metadata$Sample_ID)) {
  stop("Sample order/names in counts do not exactly match metadata$Sample_ID.")
}

metadata$Group <- factor(metadata$Group, levels = c(reference_group, comparison_group))

dge <- edgeR::DGEList(counts = counts)
keep <- edgeR::filterByExpr(dge, group = metadata$Group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ Group, data = metadata)
dge <- edgeR::estimateDisp(dge, design)
fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
test <- edgeR::glmQLFTest(fit, coef = paste0("Group", comparison_group))

edgeR_results <- edgeR::topTags(test, n = Inf)$table |>
  tibble::rownames_to_column(gene_column) |>
  dplyr::mutate(
    Regulation = dplyr::case_when(
      FDR < pvalue_threshold & logFC > logfc_threshold ~ "Up",
      FDR < pvalue_threshold & logFC < -logfc_threshold ~ "Down",
      TRUE ~ "Not significant"
    )
  )

dir.create("results/DEG", recursive = TRUE, showWarnings = FALSE)
write.csv(edgeR_results, "results/DEG/edgeR_Case_vs_Control.csv", row.names = FALSE)

p <- create_volcano_plot(
  edgeR_results,
  gene_column = gene_column,
  logfc_column = "logFC",
  pvalue_column = "FDR",
  logfc_threshold = logfc_threshold,
  pvalue_threshold = pvalue_threshold,
  y_max = (volcano_y_max), 
  title = "Volcano Plot: edgeR Case vs Control"
)
print(p)
ggplot2::ggsave("results/DEG/edgeR_volcano.png", p, width = 8, height = 6, dpi = 300)
print(table(edgeR_results$Regulation))

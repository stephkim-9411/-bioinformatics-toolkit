# ============================================================
# 02C. Compare edgeR and DESeq2 Results
# Compares fold-change estimates and identifies consensus DEGs.
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

pvalue_threshold <- 0.05
logfc_threshold <- 1

edgeR_file <- "results/DEG/edgeR_Case_vs_Control.csv"
DESeq2_file <- "results/DEG/DESeq2_Case_vs_Control.csv"

# ------------------------------------------------------------
# 2. Read edgeR and DESeq2 results
# ------------------------------------------------------------

edgeR_results <- data.table::fread(edgeR_file)

DESeq2_results <- data.table::fread(DESeq2_file)

# ------------------------------------------------------------
# 3. Combine edgeR and DESeq2 results by gene ID
# ------------------------------------------------------------

comparison <- edgeR_results |>
  dplyr::select(
    Gene_ID,
    edgeR_logFC = logFC,
    edgeR_FDR = FDR
  ) |>
  dplyr::inner_join(
    DESeq2_results |>
      dplyr::select(
        Gene_ID,
        DESeq2_logFC = log2FoldChange,
        DESeq2_FDR = padj
      ),
    by = "Gene_ID"
  )

# ------------------------------------------------------------
# 4. Identify significant and consensus DEGs
# ------------------------------------------------------------

comparison <- comparison |>
  dplyr::mutate(
    edgeR_significant =
      !is.na(edgeR_FDR) &
      edgeR_FDR < pvalue_threshold &
      abs(edgeR_logFC) > logfc_threshold,
    
    DESeq2_significant =
      !is.na(DESeq2_FDR) &
      DESeq2_FDR < pvalue_threshold &
      abs(DESeq2_logFC) > logfc_threshold,
    
    Consensus =
      edgeR_significant &
      DESeq2_significant
  )

# Select genes identified as significant by both methods
consensus_DEGs <- comparison |>
  dplyr::filter(Consensus)

message("Genes compared: ", nrow(comparison))
message("Consensus DEGs: ", nrow(consensus_DEGs))

# ------------------------------------------------------------
# 5. Save comparison and consensus DEG results
# ------------------------------------------------------------

dir.create(
  "results/DEG",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  comparison,
  "results/DEG/edgeR_vs_DESeq2_comparison.csv",
  row.names = FALSE
)

write.csv(
  consensus_DEGs,
  "results/DEG/consensus_DEGs.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Compare fold-change estimates
# ------------------------------------------------------------

p_comparison <- ggplot2::ggplot(
  comparison,
  ggplot2::aes(
    x = edgeR_logFC,
    y = DESeq2_logFC
  )
) +
  ggplot2::geom_point(
    alpha = 0.6
  ) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "edgeR vs DESeq2 Fold-Change Estimates",
    x = "edgeR log2 Fold Change",
    y = "DESeq2 log2 Fold Change"
  )

print(p_comparison)

# ------------------------------------------------------------
# 7. Save the fold-change comparison plot
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = "results/DEG/edgeR_vs_DESeq2_logFC.png",
  plot = p_comparison,
  width = 6,
  height = 6,
  dpi = 300
)

message("edgeR and DESeq2 comparison completed.")
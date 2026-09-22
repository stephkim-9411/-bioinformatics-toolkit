# ============================================================
# 00. Install Required Packages
# ============================================================

cran_packages <- c("data.table", "dplyr", "ggplot2", "ggrepel", "pheatmap")
bioc_packages <- c(
  "edgeR", "DESeq2", "clusterProfiler", "enrichplot",
  "STRINGdb", "AnnotationDbi",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"
)

cran_missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(cran_missing) > 0) install.packages(cran_missing)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_missing <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(bioc_missing) > 0) BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)

message("Package check complete.")

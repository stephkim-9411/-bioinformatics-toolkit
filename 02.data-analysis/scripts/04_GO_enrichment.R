
# ============================================================
# 04. GO Over-Representation Analysis
# BP, MF, and CC
# ============================================================

library(data.table)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

organism <- "mouse"
p_adjust_threshold <- 0.05
logfc_threshold <- 1
input_file <- "results/annotation/annotated_DEG.csv"

# ------------------------------------------------------------
# 2. Select organism database
# ------------------------------------------------------------

if (organism == "human") {
  organism_db <- org.Hs.eg.db::org.Hs.eg.db
} else if (organism == "mouse") {
  organism_db <- org.Mm.eg.db::org.Mm.eg.db
} else if (organism == "rat") {
  organism_db <- org.Rn.eg.db::org.Rn.eg.db
} else {
  stop("Unsupported organism. Use 'human', 'mouse', or 'rat'.")
}

# ------------------------------------------------------------
# 3. Read annotated DEG results
# ------------------------------------------------------------

x <- data.table::fread(input_file)

# ------------------------------------------------------------
# 4. Select significant genes
# ------------------------------------------------------------

sig <- x |>
  dplyr::filter(
    !is.na(padj),
    padj < p_adjust_threshold,
    abs(log2FoldChange) > logfc_threshold
  )

genes <- unique(
  as.character(stats::na.omit(sig$ENTREZID))
)

if (length(genes) == 0) {
  stop("No significant annotated genes available for GO enrichment.")
}

message("Number of significant genes used: ", length(genes))

# ------------------------------------------------------------
# 5. Create output directory
# ------------------------------------------------------------

dir.create(
  "results/GO",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 6A. GO Biological Process (BP)
# ============================================================

go_bp <- clusterProfiler::enrichGO(
  gene          = genes,
  OrgDb         = organism_db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

go_bp_results <- as.data.frame(go_bp)

write.csv(
  go_bp_results,
  "results/GO/GO_BP_results.csv",
  row.names = FALSE
)

if (nrow(go_bp_results) > 0) {
  p_bp <- enrichplot::dotplot(
    go_bp,
    showCategory = 15
  ) +
    ggplot2::ggtitle("GO Biological Process")
  
  print(p_bp)
  
  ggplot2::ggsave(
    filename = "results/GO/GO_BP_dotplot.png",
    plot = p_bp,
    width = 8,
    height = 6,
    dpi = 300
  )
} else {
  message("No significant GO Biological Process terms were identified.")
}

# ============================================================
# 6B. GO Molecular Function (MF)
# ============================================================

go_mf <- clusterProfiler::enrichGO(
  gene          = genes,
  OrgDb         = organism_db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

go_mf_results <- as.data.frame(go_mf)

write.csv(
  go_mf_results,
  "results/GO/GO_MF_results.csv",
  row.names = FALSE
)

if (nrow(go_mf_results) > 0) {
  p_mf <- enrichplot::dotplot(
    go_mf,
    showCategory = 15
  ) +
    ggplot2::ggtitle("GO Molecular Function")
  
  print(p_mf)
  
  ggplot2::ggsave(
    filename = "results/GO/GO_MF_dotplot.png",
    plot = p_mf,
    width = 8,
    height = 8,
    dpi = 300
  )
} else {
  message("No significant GO Molecular Function terms were identified.")
}

# ============================================================
# 6C. GO Cellular Component (CC)
# ============================================================

go_cc <- clusterProfiler::enrichGO(
  gene          = genes,
  OrgDb         = organism_db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

go_cc_results <- as.data.frame(go_cc)

write.csv(
  go_cc_results,
  "results/GO/GO_CC_results.csv",
  row.names = FALSE
)

if (nrow(go_cc_results) > 0) {
  p_cc <- enrichplot::dotplot(
    go_cc,
    showCategory = 15
  ) +
    ggplot2::ggtitle("GO Cellular Component")
  
  print(p_cc)
  
  ggplot2::ggsave(
    filename = "results/GO/GO_CC_dotplot.png",
    plot = p_cc,
    width = 8,
    height = 6,
    dpi = 300
  )
} else {
  message("No significant GO Cellular Component terms were identified.")
}


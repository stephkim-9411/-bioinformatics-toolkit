# ============================================================
# 07. Gene Set Enrichment Analysis (GSEA)
# Uses ALL ranked genes, not only significant DEGs.
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
# 3. Read all genes
# ------------------------------------------------------------

x <- data.table::fread(input_file) |>
  dplyr::filter(
    !is.na(ENTREZID),
    !is.na(log2FoldChange)
  )

# ------------------------------------------------------------
# 4. Remove duplicated Entrez IDs
# Keep the gene with the largest absolute fold change
# ------------------------------------------------------------

rank_df <- x |>
  dplyr::arrange(
    dplyr::desc(abs(log2FoldChange))
  ) |>
  dplyr::distinct(
    ENTREZID,
    .keep_all = TRUE
  )

# ------------------------------------------------------------
# 5. Create ranked gene list
# ------------------------------------------------------------

gene_list <- rank_df$log2FoldChange

names(gene_list) <- as.character(rank_df$ENTREZID)

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

# Create output directory
dir.create(
  "results/GSEA",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 6A. GSEA: GO Biological Process (BP)
# ============================================================

gsea_bp <- clusterProfiler::gseGO(
  geneList      = gene_list,
  OrgDb         = organism_db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

gsea_bp_results <- as.data.frame(gsea_bp)

write.csv(
  gsea_bp_results,
  "results/GSEA/GSEA_GO_BP_results.csv",
  row.names = FALSE
)

if (nrow(gsea_bp_results) > 0) {
  p_bp <- enrichplot::dotplot(
    gsea_bp,
    showCategory = 15,
    split = ".sign"
  ) +
    ggplot2::facet_grid(. ~ .sign) +
    ggplot2::ggtitle("GSEA: GO Biological Process")
  
  print(p_bp)
  
  ggplot2::ggsave(
    filename = "results/GSEA/GSEA_GO_BP_dotplot.png",
    plot = p_bp,
    width = 9,
    height = 7,
    dpi = 300
  )
} else {
  message("No significant GSEA GO Biological Process terms were identified.")
}

# ============================================================
# 6B. GSEA: GO Molecular Function (MF)
# ============================================================

gsea_mf <- clusterProfiler::gseGO(
  geneList      = gene_list,
  OrgDb         = organism_db,
  keyType       = "ENTREZID",
  ont           = "MF",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

gsea_mf_results <- as.data.frame(gsea_mf)

write.csv(
  gsea_mf_results,
  "results/GSEA/GSEA_GO_MF_results.csv",
  row.names = FALSE
)

if (nrow(gsea_mf_results) > 0) {
  p_mf <- enrichplot::dotplot(
    gsea_mf,
    showCategory = 15,
    split = ".sign"
  ) +
    ggplot2::facet_grid(. ~ .sign) +
    ggplot2::ggtitle("GSEA: GO Molecular Function")
  
  print(p_mf)
  
  ggplot2::ggsave(
    filename = "results/GSEA/GSEA_GO_MF_dotplot.png",
    plot = p_mf,
    width = 9,
    height = 7,
    dpi = 300
  )
} else {
  message("No significant GSEA GO Molecular Function terms were identified.")
}

# ============================================================
# 6C. GSEA: GO Cellular Component (CC)
# ============================================================

gsea_cc <- clusterProfiler::gseGO(
  geneList      = gene_list,
  OrgDb         = organism_db,
  keyType       = "ENTREZID",
  ont           = "CC",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

gsea_cc_results <- as.data.frame(gsea_cc)

write.csv(
  gsea_cc_results,
  "results/GSEA/GSEA_GO_CC_results.csv",
  row.names = FALSE
)

if (nrow(gsea_cc_results) > 0) {
  p_cc <- enrichplot::dotplot(
    gsea_cc,
    showCategory = 15,
    split = ".sign"
  ) +
    ggplot2::facet_grid(. ~ .sign) +
    ggplot2::ggtitle("GSEA: GO Cellular Component")
  
  print(p_cc)
  
  ggplot2::ggsave(
    filename = "results/GSEA/GSEA_GO_CC_dotplot.png",
    plot = p_cc,
    width = 9,
    height = 7,
    dpi = 300
  )
} else {
  message("No significant GSEA GO Cellular Component terms were identified.")
}

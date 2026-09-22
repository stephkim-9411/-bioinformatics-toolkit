# ============================================================
# 05. KEGG Over-Representation Analysis
# Identifies KEGG pathways enriched among significant DEGs.
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
# 2. Select KEGG organism code
# ------------------------------------------------------------

kegg_code <- switch(
  organism,
  human = "hsa",
  mouse = "mmu",
  rat   = "rno",
  stop("Unsupported organism. Use 'human', 'mouse', or 'rat'.")
)

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

# Extract unique Entrez gene IDs
genes <- unique(
  as.character(stats::na.omit(sig$ENTREZID))
)

if (length(genes) == 0) {
  stop("No significant annotated genes available for KEGG enrichment.")
}

message("Number of significant genes used: ", length(genes))

# ------------------------------------------------------------
# 5. Perform KEGG pathway enrichment
# ------------------------------------------------------------

kegg_results <- clusterProfiler::enrichKEGG(
  gene          = genes,
  organism      = kegg_code,
  keyType       = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20
)

# Convert the enrichment result to a data frame
kegg_results_df <- as.data.frame(kegg_results)

# ------------------------------------------------------------
# 6. Create output directory and save results
# ------------------------------------------------------------

dir.create(
  "results/KEGG",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  kegg_results_df,
  "results/KEGG/KEGG_results.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Create and save KEGG enrichment dot plot
# ------------------------------------------------------------

if (nrow(kegg_results_df) > 0) {
  p_kegg <- enrichplot::dotplot(
    kegg_results,
    showCategory = 15
  ) +
    ggplot2::ggtitle("KEGG Pathway Enrichment")
  
  print(p_kegg)
  
  ggplot2::ggsave(
    filename = "results/KEGG/KEGG_dotplot.png",
    plot = p_kegg,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message("KEGG enrichment analysis completed.")
} else {
  message("No significant KEGG pathways were identified.")
}


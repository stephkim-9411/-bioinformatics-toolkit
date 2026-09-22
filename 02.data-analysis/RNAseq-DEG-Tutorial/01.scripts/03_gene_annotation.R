# ============================================================
# 03. Gene Annotation
# ============================================================
library(data.table)
library(dplyr)
library(clusterProfiler)


# USER SETTINGS
organism <- "mouse"       # "human", "mouse", or "rat"
gene_id_type <- "SYMBOL"  # e.g. "SYMBOL", "ENSEMBL", "ENTREZID"
input_file <- "results/DEG/DESeq2_Case_vs_Control.csv"
gene_column <- "Gene_ID"

if (organism == "human") {
  organism_db <- org.Hs.eg.db::org.Hs.eg.db
  kegg_code <- "hsa"
  string_species <- 9606
} else if (organism == "mouse") {
  organism_db <- org.Mm.eg.db::org.Mm.eg.db
  kegg_code <- "mmu"
  string_species <- 10090
} else if (organism == "rat") {
  organism_db <- org.Rn.eg.db::org.Rn.eg.db
  kegg_code <- "rno"
  string_species <- 10116
} else {
  stop("organism must be 'human', 'mouse', or 'rat'.")
}

deg <- data.table::fread(input_file)

annotation <- clusterProfiler::bitr(
  unique(deg[[gene_column]]),
  fromType = gene_id_type,
  toType = unique(c("SYMBOL", "ENTREZID")),
  OrgDb = organism_db
)

names(annotation)[names(annotation) == gene_id_type] <- gene_column
annotated <- dplyr::left_join(deg, annotation, by = gene_column)

dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
write.csv(annotated, "results/annotation/annotated_DEG.csv", row.names = FALSE)

settings <- data.frame(
  organism = organism, kegg_code = kegg_code, string_species = string_species,
  gene_id_type = gene_id_type
)
write.csv(settings, "results/annotation/organism_settings.csv", row.names = FALSE)
message("Annotation complete. Original gene IDs were retained.")

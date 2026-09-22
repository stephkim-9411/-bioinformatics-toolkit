# ============================================================
# 06. STRING Protein-Protein Interaction Network
# Identifies interactions among significant DEGs.
# ============================================================

library(data.table)
library(dplyr)
library(STRINGdb)

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

organism <- "mouse"
p_adjust_threshold <- 0.05
logfc_threshold <- 1
input_file <- "results/annotation/annotated_DEG.csv"

# ------------------------------------------------------------
# 2. Select STRING species ID
# ------------------------------------------------------------

string_species <- switch(
  organism,
  human = 9606,
  mouse = 10090,
  rat   = 10116,
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
    abs(log2FoldChange) > logfc_threshold,
    !is.na(Gene_ID)
  ) |>
  dplyr::select(
    SYMBOL = Gene_ID,
    log2FoldChange,
    padj
  )

if (nrow(sig) == 0) {
  stop("No significant genes available for STRING analysis.")
}

message("Number of significant genes: ", nrow(sig))

# ------------------------------------------------------------
# 5. Connect to the STRING database
# A score threshold of 400 represents medium confidence.
# ------------------------------------------------------------

string_db <- STRINGdb::STRINGdb$new(
  species = string_species,
  score_threshold = 400
)

# ------------------------------------------------------------
# 6. Map gene symbols to STRING protein IDs
# ------------------------------------------------------------

mapped <- string_db$map(
  sig,
  "SYMBOL",
  removeUnmappedRows = TRUE
)

if (nrow(mapped) == 0) {
  stop("None of the significant genes could be mapped to STRING IDs.")
}

message("Number of genes mapped to STRING: ", nrow(mapped))

# ------------------------------------------------------------
# 7. Save mapped gene results
# ------------------------------------------------------------

dir.create(
  "results/STRING",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  mapped,
  "results/STRING/STRING_mapped_genes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Create and save the STRING interaction network
# ------------------------------------------------------------

png(
  filename = "results/STRING/STRING_network.png",
  width = 1800,
  height = 1400,
  res = 200
)

string_db$plot_network(mapped$STRING_id)

title(
  main = "STRING Protein–Protein Interaction Network"
)

dev.off()

message("STRING analysis completed.")

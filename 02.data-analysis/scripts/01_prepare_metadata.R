# ============================================================
# 01. Prepare Metadata from Sample Column Names
# ============================================================
# 01. Prepare Sample Metadata
# Creates a metadata file from sample names in the count matrix.
#
# Expected sample names:
# Control_1, Control_2, Case_1, Case_2
# ============================================================

library(data.table)
library(dplyr)

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

counts_file <- "data/example_counts.csv"
metadata_file <- "data/example_metadata.csv"
gene_column <- "Gene_ID"

reference_group <- "Control"
comparison_group <- "Case"

# ------------------------------------------------------------
# 2. Read the count matrix
# ------------------------------------------------------------

counts_raw <- data.table::fread(counts_file)

# ------------------------------------------------------------
# 3. Confirm that the gene ID column exists
# ------------------------------------------------------------

if (!gene_column %in% names(counts_raw)) {
  stop(
    paste0(
      "The gene column '",
      gene_column,
      "' was not found in the count file."
    )
  )
}

# ------------------------------------------------------------
# 4. Extract sample names
# All columns except Gene_ID are treated as samples.
# ------------------------------------------------------------

sample_names <- setdiff(
  colnames(counts_raw),
  gene_column
)

if (length(sample_names) == 0) {
  stop("No sample columns were found in the count file.")
}

# ------------------------------------------------------------
# 5. Create the sample metadata
# Sample names beginning with Control are assigned to Control.
# Sample names beginning with Case are assigned to Case.
# ------------------------------------------------------------

sample_info <- data.frame(
  Sample_ID = sample_names
) |>
  dplyr::mutate(
    Group = dplyr::case_when(
      grepl(
        "^Control",
        Sample_ID,
        ignore.case = TRUE
      ) ~ reference_group,
      
      grepl(
        "^Case",
        Sample_ID,
        ignore.case = TRUE
      ) ~ comparison_group,
      
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# 6. Check whether every sample was assigned to a group
# ------------------------------------------------------------

if (anyNA(sample_info$Group)) {
  unclassified_samples <- sample_info |>
    dplyr::filter(is.na(Group)) |>
    dplyr::pull(Sample_ID)
  
  stop(
    paste0(
      "The following samples could not be classified: ",
      paste(unclassified_samples, collapse = ", "),
      ". Edit the naming rules in 01_prepare_metadata.R."
    )
  )
}

# ------------------------------------------------------------
# 7. Set the reference and comparison group order
# ------------------------------------------------------------

sample_info$Group <- factor(
  sample_info$Group,
  levels = c(
    reference_group,
    comparison_group
  )
)

# ------------------------------------------------------------
# 8. Create the output directory
# ------------------------------------------------------------

dir.create(
  "data",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 9. Save the metadata file
# ------------------------------------------------------------

write.csv(
  sample_info,
  file = metadata_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Display the completed metadata
# ------------------------------------------------------------

print(sample_info)

message(
  "Sample metadata successfully saved to: ",
  metadata_file
)
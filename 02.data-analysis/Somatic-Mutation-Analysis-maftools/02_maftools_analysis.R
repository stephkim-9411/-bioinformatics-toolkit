# maftools: mutation and optional CNV analysis
# Input preparation: 01_MAF_preparation_and_annotation.md
# Run from the project root. Edit settings below, then run top to bottom.
# Default input: public TCGA AML example installed with maftools.
# One-time package installation is documented in the preparation guide.

library(maftools)
library(data.table)
library(dplyr)
library(ggplot2)


# 4. Set paths and optional analyses

# Save this script at your project root, with input files inside data/. Set your working directory to that folder before running the script.

use_example <- TRUE
maf_file <- "data/cohort.maf"  # Or data/oncokb_annotated.maf
output_dir <- "results/maftools_example" # Change for each dataset/run

use_cnv <- FALSE
cnv_file <- "data/cnv.tsv"
use_oncokb <- FALSE           # TRUE only for an already annotated MAF
reference_build <- "hg19"    # Bundled AML example; use hg38 for GRCh38 data

top_genes <- 20L
exclude_samples <- character(0) # For example: c("A16")
selected_genes <- c("NPM1", "FLT3", "DNMT3A", "IDH1", "IDH2", "RUNX1", "TP53")
run_rainfall <- TRUE
run_signatures <- FALSE      # Advanced and computationally intensive
chosen_signature_n <- NULL   # Choose after inspecting rank diagnostics

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (use_example) {
  maf_file <- "data/tcga_laml.maf.gz"
  if (!file.exists(maf_file)) {
    maf_file <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
  }
  if (use_cnv || use_oncokb) stop("Set use_example = FALSE for your CNV/OncoKB files.")
  reference_build <- "hg19"
}
stopifnot(file.exists(maf_file), reference_build %in% c("hg19", "hg38"))

# 5. Read and check the MAF

# `read.maf()` builds a **MAF object**, which contains processed variants and summaries. Its `@data` slot holds variants classified as nonsynonymous by maftools; other classes are in `@maf.silent`. `useAll = TRUE` means that `Mutation_Status` is not used to filter rows; it does **not** mean all consequences appear in nonsynonymous plots. Only supply the intended, quality-filtered somatic dataset.

maf_table <- data.table::fread(maf_file, skip = "Hugo_Symbol", sep = "\t")
# Normalize the maftools AML example's column names without changing its variants.
if ("End_position" %in% names(maf_table) && !"End_Position" %in% names(maf_table)) {
  data.table::setnames(maf_table, "End_position", "End_Position")
}
if (use_example && "Protein_Change" %in% names(maf_table) &&
    !"HGVSp_Short" %in% names(maf_table)) {
  data.table::setnames(maf_table, "Protein_Change", "HGVSp_Short")
}
required <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
              "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
              "Variant_Type", "Tumor_Sample_Barcode")
missing <- setdiff(required, names(maf_table))
if (length(missing)) stop("Missing MAF columns: ", paste(missing, collapse = ", "))
if (anyNA(maf_table[, ..required]) ||
    any(vapply(maf_table[, ..required], function(x) any(trimws(as.character(x)) == ""), logical(1)))) {
  stop("Required MAF fields contain missing/blank values; inspect the input.")
}
maf_table <- maf_table %>%
  dplyr::filter(!Tumor_Sample_Barcode %in% exclude_samples)
stopifnot(nrow(maf_table) > 0)

mutation_maf <- maftools::read.maf(
  maf = maf_table, useAll = TRUE,
  removeDuplicatedVariants = TRUE, verbose = TRUE
)
mutation_maf

# Repeated transcript annotations are deduplicated by maftools. Resolve contradictory annotations upstream; keeping the first duplicate does not adjudicate its biological interpretation. Check sample IDs and genome build before continuing.

# 6. Optional: include CNV as an alteration category

# Here “CNV as a factor” means **Amp/Del categories in a combined oncoplot**. It does not mean fitting a statistical model with CNV as a covariate. The mutation-only object remains available for protein and positional analyses.

# Provide a gene-level TSV with these input columns:

# | Hugo_Symbol | Tumor_Sample_Barcode | Copy_Number_Alteration |
# |---|---|---|
# | ERBB2 | Tumor01 | Amplification |
# | CDKN2A | Tumor02 | DeepDeletion |

# These rows illustrate the format only. Segment-level CNVkit output needs gene assignment and justified copy-number calling first. Do not apply universal absolute-copy-number thresholds without considering caller definitions, purity, and ploidy.

# The mapping below retains explicit amplification/deep-deletion calls. `Gain`, shallow loss, generic `Deletion`, and unknown labels are exported for review rather than silently recoded. Add mappings only after checking what your caller means. For thresholded GISTIC data, +2/-2 can represent amplification/deep deletion; do not apply that interpretation to raw copy numbers or log2 ratios.

cn_table <- NULL
combined_maf <- NULL
if (use_cnv) {
  stopifnot(file.exists(cnv_file))
  cnv <- data.table::fread(cnv_file)
  cnv_required <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Copy_Number_Alteration")
  if (!all(cnv_required %in% names(cnv))) stop("Check the CNV column names.")
  cnv <- cnv %>%
    dplyr::filter(!Tumor_Sample_Barcode %in% exclude_samples) %>%
    dplyr::mutate(CN = dplyr::case_when(
      Copy_Number_Alteration %in% c("Amplification", "Amp") ~ "Amp",
      Copy_Number_Alteration %in% c("DeepDeletion", "Deep Deletion", "Del") ~ "Del",
      TRUE ~ NA_character_
    ))
  data.table::fwrite(dplyr::filter(cnv, is.na(CN)),
                    file.path(output_dir, "CNV_unmapped_rows.csv"))
  # Column ORDER is important: gene, sample, CN status.
  cn_table <- cnv %>%
    dplyr::filter(!is.na(CN)) %>%
    dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, CN) %>%
    dplyr::distinct()
  if (!nrow(cn_table)) stop("No retained CNVs; review labels or set use_cnv = FALSE.")
  if (anyNA(cn_table) || any(trimws(cn_table$Hugo_Symbol) == "") ||
      any(trimws(cn_table$Tumor_Sample_Barcode) == "")) stop("Missing CNV gene/sample IDs.")
  conflicts <- cn_table %>%
    dplyr::count(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    dplyr::filter(n > 1)
  if (nrow(conflicts)) stop("Conflicting Amp/Del calls for a gene/sample; resolve first.")
  cnv_only_ids <- setdiff(cn_table$Tumor_Sample_Barcode, maf_table$Tumor_Sample_Barcode)
  if (length(cnv_only_ids)) {
    stop("CNV-only IDs need review against the cohort/sample manifest: ",
         paste(cnv_only_ids, collapse = ", "))
  }
  print(table(cn_table$CN))
  print(setdiff(unique(maf_table$Tumor_Sample_Barcode), unique(cn_table$Tumor_Sample_Barcode)))
  data.table::fwrite(cn_table, file.path(output_dir, "CNV_maftools_input.tsv"), sep = "\t")
  combined_maf <- maftools::read.maf(
    maf = maf_table, cnTable = as.data.frame(cn_table),
    useAll = TRUE, removeDuplicatedVariants = TRUE
  )
}
plot_maf <- if (use_cnv) combined_maf else mutation_maf

# A sample absent from the retained CNV table may have no qualifying calls **or** may not have been assayed. Verify this with sample metadata. This simple tutorial stops on CNV-only sample IDs instead of silently changing the cohort denominator. Samples with no retained mutation or CNV require an explicit full-cohort manifest and careful denominator handling; do not add fake mutation rows.

# 7. Oncoplots: genes by samples

# Rows are genes; columns are samples; colors show alteration types. With CNV enabled, the plot also includes Amp/Del. Frequencies describe the represented MAF cohort and the selected alterations, not automatically every enrolled patient. `removeNonMutated = FALSE` retains represented samples without changes in the displayed genes.

# A small helper displays a plot in R and saves it as a PDF. The saved output is redrawn, so plotting functions should not modify data.

show_and_save <- function(filename, draw, width = 12, height = 7) {
  draw()
  grDevices::pdf(file.path(output_dir, filename), width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}



show_and_save("Oncoplot.pdf", function() {
  maftools::oncoplot(
    maf = plot_maf, top = top_genes, removeNonMutated = FALSE,
    showTumorSampleBarcodes = TRUE, SampleNamefontSize = 0.6,
    fontSize = 0.8, draw_titv = FALSE,
    titleText = if (use_cnv) "Sequence variants and copy-number alterations" else "Sequence variants"
  )
})

# Optional gene panel and sample order

# Change `selected_genes` in the settings. Samples can also be manually ordered with `sample_order <- c("A1", "A2", "A10")`; include all intended represented samples. Alphabetical order below is deliberately simple.

plot_genes <- intersect(selected_genes, as.character(maftools::getGeneSummary(plot_maf)$Hugo_Symbol))
sample_order <- sort(as.character(maftools::getSampleSummary(plot_maf)$Tumor_Sample_Barcode))
if (length(plot_genes)) {
  show_and_save("Selected_genes_oncoplot.pdf", function() {
    maftools::oncoplot(maf = plot_maf, genes = plot_genes,
      sampleOrder = sample_order, keepGeneOrder = TRUE,
      removeNonMutated = FALSE, showTumorSampleBarcodes = TRUE,
      SampleNamefontSize = 0.6, draw_titv = FALSE)
  })
}

# 8. Mutation summary dashboard and tables

# The dashboard summarizes mutation consequences, types, and per-sample counts. It uses the mutation-only object so CNV events are not mixed into sequence mutation counts. Counts are not tumor mutational burden (TMB): TMB requires an appropriate callable megabase denominator and filtering definition. `rmOutlier = TRUE` affects the display, not the underlying input.

show_and_save("Mutation_MAF_summary.pdf", function() {
  maftools::plotmafSummary(maf = mutation_maf, rmOutlier = TRUE,
    addStat = "median", dashboard = TRUE, titvRaw = FALSE)
})
gene_summary <- maftools::getGeneSummary(mutation_maf)
sample_summary <- maftools::getSampleSummary(mutation_maf)
data.table::fwrite(gene_summary, file.path(output_dir, "Mutation_gene_summary.csv"))
data.table::fwrite(sample_summary, file.path(output_dir, "Mutation_sample_summary.csv"))
head(gene_summary)
head(sample_summary)
if (use_cnv) {
  data.table::fwrite(maftools::getGeneSummary(combined_maf),
                    file.path(output_dir, "Combined_gene_summary.csv"))
  data.table::fwrite(maftools::getSampleSummary(combined_maf),
                    file.path(output_dir, "Combined_sample_summary.csv"))
}

# 9. Protein lollipop plots for the top 20 genes

# Lollipop plots place mutations along a protein and can reveal recurrent positions. They need interpretable protein-change annotations and compatible protein/transcript information. CNVs do not have amino-acid positions. Genes that cannot be plotted are recorded with the error message instead of being silently skipped.

plot_log <- data.frame(Analysis = character(), Item = character(), Message = character())
if ("HGVSp_Short" %in% names(mutation_maf@data)) {
  local({
    grDevices::pdf(file.path(output_dir, "Top20_LollipopPlots.pdf"), width = 10, height = 5)
    on.exit(grDevices::dev.off())
    for (gene_name in head(gene_summary$Hugo_Symbol, top_genes)) {
      tryCatch(maftools::lollipopPlot(maf = mutation_maf, gene = gene_name,
        AACol = "HGVSp_Short", showMutationRate = TRUE), error = function(e) {
          plot_log <<- rbind(plot_log, data.frame(Analysis = "Lollipop",
            Item = gene_name, Message = conditionMessage(e)))
          message(gene_name, ": ", conditionMessage(e))
        })
    }
  })
} else {
  message("Lollipops skipped: HGVSp_Short is unavailable.")
}

# 10. Optional: OncoKB annotation and interpretation

# For upstream annotation, see 01_MAF_preparation_and_annotation.md.

# `ONCOGENIC` describes oncogenicity; `MUTATION_EFFECT` describes functional effect; `HIGHEST_LEVEL` records the highest returned therapeutic evidence level. Oncogenicity and treatment evidence are different concepts. R1/R2 indicate resistance, not sensitivity. Blank evidence is labeled **No reported level**, not “non-actionable.” Interpret annotations in their tumor context and record the annotation date/version.

# These summaries count retained **variant–sample records**, not unique patients or globally unique alleles. They intentionally cover the nonsynonymous sequence variants in `mutation_maf@data`, matching the oncoplot. They do not summarize CNV treatment evidence.

oncogenic_maf <- NULL
oncogenic_only <- NULL
if (use_oncokb) {
  if (!"ONCOGENIC" %in% names(mutation_maf@data)) stop("MAF lacks ONCOGENIC annotation.")
  oncogenic_only <- as.data.frame(mutation_maf@data) %>%
    dplyr::filter(ONCOGENIC %in% c("Oncogenic", "Likely Oncogenic"))
  data.table::fwrite(oncogenic_only,
    file.path(output_dir, "OncoKB_Oncogenic_Likely_Oncogenic_Variants.csv"))
  if (nrow(oncogenic_only)) {
    oncogenic_maf <- maftools::read.maf(maf = oncogenic_only, useAll = TRUE,
                                      removeDuplicatedVariants = TRUE)
    show_and_save("OncoKB_Oncoplot.pdf", function() {
      maftools::oncoplot(maf = oncogenic_maf, top = top_genes,
        showTumorSampleBarcodes = TRUE, SampleNamefontSize = 0.6,
        draw_titv = FALSE, titleText = "Oncogenic / Likely Oncogenic sequence variants")
    })
    show_and_save("OncoKB_MAF_Summary.pdf", function() {
      maftools::plotmafSummary(maf = oncogenic_maf, rmOutlier = TRUE,
                              addStat = "median", dashboard = TRUE, titvRaw = FALSE)
    })
    data.table::fwrite(maftools::getGeneSummary(oncogenic_maf),
                      file.path(output_dir, "OncoKB_gene_summary.csv"))
    data.table::fwrite(maftools::getSampleSummary(oncogenic_maf),
                      file.path(output_dir, "OncoKB_sample_summary.csv"))
  } else message("No Oncogenic/Likely Oncogenic nonsynonymous records; optional plots skipped.")
}

# The filtered MAF may omit samples with no qualifying records, changing its displayed denominator. Use a cohort manifest to calculate full-cohort prevalence. CNVs are not included in this filtered object: restricting CNVs to genes with oncogenic mutations does not establish the CNVs' oncogenicity. Annotate and filter each CNV independently if an oncogenic-CNV plot is needed.

# Evidence-level and mutation-effect barplots

# The original levels are preserved, including unexpected labels, so changed vocabularies do not become missing values. Each chart is saved as PDF and PNG, with its underlying CSV.

if (!is.null(oncogenic_only) && nrow(oncogenic_only)) {
  for (field in c("HIGHEST_LEVEL", "MUTATION_EFFECT")) {
    if (!field %in% names(oncogenic_only)) {
      message("Skipped missing annotation column: ", field)
      next
    }
    labels <- trimws(as.character(oncogenic_only[[field]]))
    labels[is.na(labels) | labels == ""] <- if (field == "HIGHEST_LEVEL") "No reported level" else "Unknown"
    summary_table <- data.frame(Category = labels) %>%
      dplyr::count(Category, name = "Variant_Sample_Count") %>%
      dplyr::arrange(dplyr::desc(Variant_Sample_Count))
    data.table::fwrite(summary_table, file.path(output_dir, paste0("OncoKB_", field, "_summary.csv")))
    p <- ggplot2::ggplot(summary_table,
      ggplot2::aes(x = reorder(Category, Variant_Sample_Count), y = Variant_Sample_Count)) +
      ggplot2::geom_col(fill = "#326A8A", width = 0.75) +
      ggplot2::geom_text(ggplot2::aes(label = Variant_Sample_Count), hjust = -0.15) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
      ggplot2::labs(title = paste("OncoKB:", field),
        subtitle = "Oncogenic / Likely Oncogenic sequence variants",
        x = NULL, y = "Variant–sample records") +
      ggplot2::theme_classic(base_size = 12)
    print(p)
    for (ext in c("pdf", "png")) {
      ggplot2::ggsave(file.path(output_dir, paste0("OncoKB_", field, ".", ext)),
                     p, width = 9, height = 5, dpi = 300)
    }
  }
}

# Level 1 evidence table and recurrent protein changes

# Level 1 output is a review table for records with that returned evidence level, not a patient-level treatment recommendation. Optional annotation columns are selected only when present.

if (!is.null(oncogenic_only) && nrow(oncogenic_only)) {
  if ("HIGHEST_LEVEL" %in% names(oncogenic_only)) {
    level1_table <- oncogenic_only %>%
      dplyr::filter(trimws(HIGHEST_LEVEL) == "LEVEL_1") %>%
      dplyr::select(dplyr::any_of(c("Tumor_Sample_Barcode", "Hugo_Symbol",
        "HGVSp_Short", "HGVSc", "Variant_Classification", "ONCOGENIC",
        "MUTATION_EFFECT", "HIGHEST_LEVEL", "LEVEL_1", "TX_CITATIONS"))) %>%
      dplyr::arrange(Hugo_Symbol, Tumor_Sample_Barcode)
    data.table::fwrite(level1_table, file.path(output_dir, "OncoKB_LEVEL1_variants.csv"))
  }
  if ("HGVSp_Short" %in% names(oncogenic_only)) {
    protein_summary <- oncogenic_only %>%
      dplyr::filter(!is.na(HGVSp_Short), HGVSp_Short != "") %>%
      dplyr::group_by(Hugo_Symbol, HGVSp_Short, Variant_Classification, ONCOGENIC) %>%
      dplyr::summarise(Variant_Sample_Count = dplyr::n(),
        Sample_Count = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") %>%
      dplyr::arrange(Hugo_Symbol, dplyr::desc(Sample_Count))
    data.table::fwrite(protein_summary, file.path(output_dir, "OncoKB_protein_change_summary.csv"))
  }
}

# 11. Rainfall plots: distance between mutations

# A rainfall plot shows the spacing of SNVs along the genome. Close clusters appear at low intermutation distances. Use the unfiltered mutation MAF, including its retained synonymous calls, rather than just OncoKB-selected variants. Targeted/exome coverage influences apparent spacing; these plots alone do not establish kataegis. The build controls chromosome sizes and must match the input.

if (run_rainfall) {
  local({
    grDevices::pdf(file.path(output_dir, "RainfallPlots.pdf"), width = 12, height = 6)
    on.exit(grDevices::dev.off())
    for (sample_id in unique(as.character(maf_table$Tumor_Sample_Barcode))) {
      tryCatch(maftools::rainfallPlot(maf = mutation_maf, tsb = sample_id,
        ref.build = reference_build, detectChangePoints = FALSE), error = function(e) {
          plot_log <<- rbind(plot_log, data.frame(Analysis = "Rainfall",
            Item = sample_id, Message = conditionMessage(e)))
          message(sample_id, ": ", conditionMessage(e))
        })
    }
  })
}
data.table::fwrite(plot_log, file.path(output_dir, "Plot_error_log.csv"))

# 12. Advanced optional analysis: SBS signatures

# This section builds a 96-channel substitution-context matrix and explores signature extraction. Use all suitable quality-filtered somatic SNVs, including synonymous SNVs, not an oncogenic-only subset. Restricting to driver variants distorts the mutation spectrum. Small cohorts, few mutations, or targeted panels may not support reliable de novo extraction. Similarity to a reference signature does not prove a causal repair defect.

# Install the matching reference package and NMF manually, then enable `run_signatures`. First inspect the rank diagnostics with `chosen_signature_n <- NULL`. Set a justified rank and rerun to extract signatures; four is not a universal choice.

# Optional signature dependencies: run these manually once.
# install.packages("NMF")
# # Choose ONE reference matching your data:
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# # BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")



if (run_signatures) {
  genome_package <- if (reference_build == "hg38") {
    "BSgenome.Hsapiens.UCSC.hg38"
  } else "BSgenome.Hsapiens.UCSC.hg19"
  if (!requireNamespace(genome_package, quietly = TRUE) ||
      !requireNamespace("NMF", quietly = TRUE)) stop("Install reference genome and NMF first.")
  tnm <- maftools::trinucleotideMatrix(maf = mutation_maf,
    ref_genome = genome_package, useSyn = TRUE)
  set.seed(42)
  sig_estimate <- maftools::estimateSignatures(mat = tnm,
    nMin = 2, nTry = 6, nrun = 20, parallel = 1)
  show_and_save("Signature_rank_diagnostics.pdf", function() {
    maftools::plotCophenetic(sig_estimate)
  })
  if (!is.null(chosen_signature_n)) {
    stopifnot(chosen_signature_n %in% 2:6)
    set.seed(42)
    sbs_signatures <- maftools::extractSignatures(mat = tnm,
      n = chosen_signature_n, parallel = 1)
    cosmic_comparison <- maftools::compareSignatures(nmfRes = sbs_signatures, sig_db = "SBS")
    show_and_save("SBS_signatures.pdf", function() {
      maftools::plotSignatures(nmfRes = sbs_signatures, contributions = FALSE)
    })
    show_and_save("SBS_contributions.pdf", function() {
      maftools::plotSignatures(nmfRes = sbs_signatures, contributions = TRUE)
    })
    saveRDS(list(signatures = sbs_signatures, comparison = cosmic_comparison),
            file.path(output_dir, "SBS_results.rds"))
  }
}

# 13. Save reproducibility information

saveRDS(list(mutations = mutation_maf, combined = combined_maf, oncogenic = oncogenic_maf),
        file.path(output_dir, "MAF_objects.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))



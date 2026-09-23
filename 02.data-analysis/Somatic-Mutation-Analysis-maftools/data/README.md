# Example data

`tcga_laml.maf.gz` is the public TCGA acute myeloid leukemia (AML/LAML) example distributed by maftools. It is real research data, not synthetic data. The downloaded bytes are unchanged.

- Source: https://github.com/PoisonAlien/maftools/blob/master/inst/extdata/tcga_laml.maf.gz
- Downloaded: 2026-09-23
- Genome build recorded in the file: 37 (use hg19 in this tutorial).
- Raw input: 2,207 rows and 193 distinct `Tumor_Sample_Barcode` values. These are input counts, not expected nonsynonymous output counts.
- SHA-256: `d102b071a052265b6f8ad7947bad1d58d3e3036fd17d6b274f7ea09a376cd6a0`
- The tutorial normalizes `End_position` to `End_Position` and `Protein_Change` to `HGVSp_Short` in memory.
- Upstream maftools license text is retained in `maftools_LICENSE.txt`; upstream attribution and data provenance are preserved.

This MAF does not contain CNV calls or OncoKB annotations. Do not fabricate these annotations for the public example. The optional sections require appropriately matched external inputs.

Use `citation("maftools")` in R for the package citation. See the official [maftools vignette](https://github.com/PoisonAlien/maftools/blob/master/vignettes/maftools.Rmd) for the AML data context and original study reference.

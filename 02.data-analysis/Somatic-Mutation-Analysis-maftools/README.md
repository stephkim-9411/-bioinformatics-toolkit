# maftools tutorial: MAF preparation, mutation analysis, and optional CNV

A step-by-step R Markdown tutorial using the **public TCGA AML example supplied by maftools**. No private input files are needed for the default workflow. The example is real public research data, not simulated data.

## Files

| File | Purpose |
|---|---|
| `01_maftools_from_MAF_to_CNV.Rmd` | Explanations and analysis code |
| `data/tcga_laml.maf.gz` | Unmodified maftools example MAF |
| `data/README.md` | Data source, input counts, genome build, and checksum |
| `data/maftools_LICENSE.txt` | Upstream license notice |

## Run the example

1. Download this folder and open the Rmd in RStudio.
2. Run the installation chunk once.
3. Keep `use_example <- TRUE`, `use_cnv <- FALSE`, and `use_oncokb <- FALSE`.
4. Click **Knit → Knit to HTML**. Output files are created in `results/maftools_example/`.

The tutorial reads the included example file. If it is missing, it falls back to the copy installed with maftools. A different installed package version may contain a different example file.

## Analyses

- VCF-to-MAF preparation instructions and combining per-sample MAFs.
- Top-gene and selected-gene oncoplots.
- Mutation summary dashboard and gene/sample summary tables.
- Top-gene protein lollipop plots.
- Sample-level rainfall plots.
- Optional CNV Amp/Del integration.
- Optional OncoKB evidence/effect plots and variant tables.
- Optional SBS signature analysis.

## Expected results

Default outputs include `Oncoplot.pdf`, `Selected_genes_oncoplot.pdf`, `Mutation_MAF_summary.pdf`, gene/sample CSVs, `Top20_LollipopPlots.pdf`, `RainfallPlots.pdf`, `Plot_error_log.csv`, `MAF_objects.rds`, and `sessionInfo.txt`.

The MAF contains 2,207 input rows from 193 sample IDs. maftools filters consequences and deduplicates records, so output mutation counts can differ. Lollipop and rainfall failures are logged. No pre-rendered results are included: R was unavailable in the preparation environment, so this package has not been executed or knitted there.

## Use your own data

Set `use_example <- FALSE`, update `maf_file` and `output_dir`, and set the correct genome build. Enable `use_cnv` only with matched gene-level CNV calls, and `use_oncokb` only with an OncoKB-annotated MAF. The AML example has neither CNV calls nor OncoKB annotation, so those sections remain disabled in the example run.

## Publish on GitHub

Upload the contents of this folder to your chosen tutorial directory, preserving `data/`. The included public example should remain separate from private analysis inputs. Run and inspect the tutorial locally before adding a rendered HTML or figures as expected results.

## Acknowledgments

Data and software: [maftools](https://github.com/PoisonAlien/maftools). Use `citation("maftools")` for the software citation. Data provenance and the upstream license notice are in `data/`.

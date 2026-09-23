# MAF preparation and annotation

Prepare your input here, then run `02_maftools_analysis.R`. Skip preparation if using the public TCGA AML example distributed with maftools (real research data, not simulated data).

MAF means Mutation Annotation Format. maftools analyzes MAF files; it does not call variants or run VEP or OncoKB annotation.

## 1. Start from a VCF: make a MAF

If you already have an annotated MAF, go directly to `02_maftools_analysis.R`. A VCF cannot be converted correctly by changing its extension: gene consequences and transcript selection must be handled during annotation.

One established workflow uses **vcf2maf with Ensembl VEP**. Install vcf2maf, VEP, the appropriate VEP cache, and a matching reference FASTA using their official instructions. The following terminal command is a template, run separately from the analysis script. Replace all paths and sample names. The VCF must already contain your quality-filtered somatic calls; conversion is not a somatic filtering step.

Create `data/per_sample_maf/` before running the conversion.

```bash
perl /path/to/vcf2maf/vcf2maf.pl \
  --input-vcf data/Tumor01.filtered.vcf \
  --output-maf data/per_sample_maf/Tumor01.maf \
  --tumor-id Tumor01 \
  --normal-id Normal01 \
  --ref-fasta /path/to/GRCh38.fa \
  --ncbi-build GRCh38 \
  --vep-path /path/to/ensembl-vep \
  --vep-data /path/to/vep-cache
```

Use matching VCF genotype sample names. If they are literally `TUMOR` and `NORMAL`, add `--vcf-tumor-id TUMOR --vcf-normal-id NORMAL` while retaining your desired output IDs. Omit normal-related options for a tumor-only input. Keep the VCF, FASTA, and annotation cache on the same genome assembly. Retain conversion logs and inspect reference-mismatch warnings.

Record the exact converter/VEP/cache versions used for reproducibility.

### Combine one MAF per sample

After converting each sample, run this chunk manually. Use a dedicated folder containing only the per-sample MAFs; the output goes outside that folder to avoid re-importing it on a later run.

```r
maf_files <- list.files("data/per_sample_maf", pattern = "\\.maf$", full.names = TRUE)
stopifnot(length(maf_files) > 0)
maf_list <- lapply(maf_files, function(path) {
  data.table::fread(path, skip = "Hugo_Symbol", sep = "\t")
})
cohort_maf <- data.table::rbindlist(maf_list, use.names = TRUE, fill = TRUE)
data.table::fwrite(cohort_maf, "data/cohort.maf", sep = "\t", quote = FALSE, na = "")
```

Do not invent consequence labels when converting an arbitrary spreadsheet. Obtain those annotations first. A minimal maftools input needs these fields (a complete standard MAF contains more):

| Column | Meaning |
|---|---|
| `Hugo_Symbol` | Gene symbol |
| `Chromosome` | Chromosome |
| `Start_Position`, `End_Position` | Genomic coordinates using MAF conventions |
| `Reference_Allele` | Reference allele |
| `Tumor_Seq_Allele2` | Tumor alternate allele |
| `Variant_Classification` | MAF consequence, e.g. `Missense_Mutation` |
| `Variant_Type` | SNP, INS, DEL, etc. |
| `Tumor_Sample_Barcode` | Exact sample identifier |

Keep `NCBI_Build` and `HGVSp_Short` when available. Protein plots use `HGVSp_Short`; it must match the chosen protein transcript. For the bundled AML example, the existing `Protein_Change` column is renamed to `HGVSp_Short` in memory.

## 2. Optional: add OncoKB annotations

OncoKB annotation is a separate upstream operation using the official annotator and authorized API access. Follow its README to install dependencies and obtain access. Do not commit an API token to GitHub. Supply the correct tumor type (OncoTree code) and reference assembly. Example terminal template:

```bash
python /path/to/oncokb-annotator/MafAnnotator.py \
  -i data/cohort.maf -o data/oncokb_annotated.maf \
  -b "$ONCOKB_API_TOKEN" -t YOUR_ONCOTREE_CODE -r GRCh38
```

## 3. Run the analysis script

Install the analysis packages once in R:

```r
install.packages(c("BiocManager", "data.table", "dplyr", "ggplot2", "R.utils"))
BiocManager::install("maftools")
```

Open `02_maftools_analysis.R` in RStudio and set your working directory to the project folder containing `data/`. Run the script from top to bottom.

- **Example data:** leave `use_example <- TRUE`. The script locates the AML MAF installed with maftools; no separate download is needed.
- **Your MAF:** set `use_example <- FALSE` and `maf_file <- "data/cohort.maf"`.
- **OncoKB-annotated MAF:** also set `maf_file <- "data/oncokb_annotated.maf"` and `use_oncokb <- TRUE`.
- **CNV:** set `use_cnv <- TRUE` and supply the gene-level TSV described in the script.

Results are saved under `output_dir`. Change that setting for your own dataset. The script includes oncoplots, mutation summaries, lollipop and rainfall plots, optional OncoKB summaries and CNVs, and optional SBS signatures. Review `Plot_error_log.csv` for skipped plots.

The script was separated from the supplied tutorial and checked structurally, but not executed here because R is unavailable. Annotation commands require your own installations, inputs, and authorized OncoKB access.

## References

- [maftools](https://github.com/poisonalien/maftools)
- [vcf2maf installation](https://github.com/mskcc/vcf2maf)
- [OncoKB annotator](https://github.com/oncokb/oncokb-annotator)

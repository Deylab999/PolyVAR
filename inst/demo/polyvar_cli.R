#!/usr/bin/env Rscript
# polyvar_cli.R  --  Command-line interface for PolyVAR
# Usage: polyvar --bfile <prefix> --pheno <file> --method polyvar [options]

suppressPackageStartupMessages({
  library(optparse)
  library(polyvar)
})

option_list <- list(
  make_option("--bfile",
    type    = "character", default = NULL,
    help    = "PLINK file prefix (no extension) [required]"),
  make_option("--pheno",
    type    = "character", default = NULL,
    help    = "Phenotype file: FID IID pheno1 pheno2 ... [required]"),
  make_option("--out",
    type    = "character", default = "polyvar_out",
    help    = "Output file prefix [default: polyvar_out]"),
  make_option("--method",
    type    = "character", default = "polyvar",
    help    = paste("vQTL method: polyvar | polyvar_lite | levene | bf |",
                    "bartlett | dglm | all [default: polyvar]")),
  make_option("--batch-size",
    type    = "integer",   default = 1000L,
    help    = "SNPs per batch [default: 1000]"),
  make_option("--alpha",
    type    = "double",    default = 5e-8,
    help    = "Genome-wide significance threshold [default: 5e-8]"),
  make_option("--min-maf",
    type    = "double",    default = 0.01,
    help    = "Minimum MAF filter [default: 0.01]"),
  make_option("--Q-lite",
    type    = "integer",   default = 100L,
    help    = "Bins for PolyVAR-lite [default: 100]"),
  make_option("--thread",
    type    = "integer",   default = 1L,
    help    = "Number of threads [default: 1, OpenMP parallelism]"),
  make_option("--verbose",
    action  = "store_true", default = TRUE,
    help    = "Print progress messages [default: TRUE]"),
  make_option("--version",
    action  = "store_true", default = FALSE,
    help    = "Print version and exit")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = paste(
    "\nPolyVAR v0.5.0 -- Contrastive Variance QTL Detection",
    "Eliminates mean-effect contamination via score-space projection.",
    "Companion paper: PolyVAR: From GWAS Basics to Contrastive vQTL Testing",
    sep = "\n"
  ),
  epilogue = paste(
    "Examples:",
    "  polyvar --bfile ukbb_chr1 --pheno nmr.tsv --method polyvar",
    "  polyvar --bfile ukbb_chr1 --pheno nmr.tsv --method all --out results",
    "  polyvar --bfile ukbb_chr1 --pheno prot.tsv --method polyvar_lite --Q-lite 50",
    sep = "\n"
  )
))

if (opt$version) {
  cat("PolyVAR version 0.5.0\n")
  quit(status = 0)
}

if (is.null(opt$bfile) || is.null(opt$pheno)) {
  cat("ERROR: --bfile and --pheno are required.\n\n")
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

# Check files exist
for (ext in c(".bed",".bim",".fam")) {
  f <- paste0(opt$bfile, ext)
  if (!file.exists(f)) stop("File not found: ", f)
}
if (!file.exists(opt$pheno)) stop("Phenotype file not found: ", opt$pheno)

methods_to_run <- if (opt$method == "all")
  c("polyvar","polyvar_lite","levene","bf","bartlett","dglm") else opt$method

cat(sprintf("\n========================================\n"))
cat(sprintf("  PolyVAR v0.5.0\n"))
cat(sprintf("  bfile  : %s\n", opt$bfile))
cat(sprintf("  pheno  : %s\n", opt$pheno))
cat(sprintf("  method : %s\n", paste(methods_to_run, collapse=", ")))
cat(sprintf("  output : %s\n", opt$out))
cat(sprintf("========================================\n\n"))

t0 <- proc.time()

for (m in methods_to_run) {
  cat(sprintf("--- Running %s ---\n", toupper(m)))
  polyvar_scan(
    bfile         = opt$bfile,
    pheno_file    = opt$pheno,
    output_prefix = opt$out,
    method        = m,
    batch_size    = opt$`batch-size`,
    alpha_thresh  = opt$alpha,
    Q_lite        = opt$`Q-lite`,
    min_maf       = opt$`min-maf`,
    verbose       = opt$verbose
  )
}

t_total <- proc.time() - t0
cat(sprintf("\nTotal wall time: %.1f seconds\n", t_total[3]))
cat("Done.\n")

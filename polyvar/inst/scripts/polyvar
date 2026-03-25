#!/usr/bin/env Rscript
# =============================================================================
#  polyvar  --  Contrastive Variance QTL Detection
#  Command-line interface, modelled on OSCA (Yang et al.)
#
#  Usage examples:
#    polyvar --bfile ukbb_chr1 --pheno nmr.phen --out results
#    polyvar --bfile ukbb_chr1 --pheno prot.phen --method polyvar_lite --out results
#    polyvar --bfile ukbb_chr1 --pheno nmr.phen --method all --out results
#    polyvar --bfile ukbb_chr1 --pheno prot.phen --method levene --chr 1 --out lev
#    polyvar --bfile ukbb_chr1 --pheno nmr.phen --extract snplist.txt --out filtered
#
#  Phenotype file format (space/tab separated, header required):
#    FID  IID  pheno1  pheno2  ...
#
#  Output: one .vqtl file per phenotype per method
#    SNP  CHR  BP  A1  A2  MAF  BETA  T_STAT  PVAL  GAMMA2
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(polyvar)
})

VERSION <- "0.5.0"

# ─────────────────────────────────────────────────────────────────────────────
# Banner
# ─────────────────────────────────────────────────────────────────────────────
print_banner <- function() {
  cat("*********************************************\n")
  cat("* PolyVAR -- Contrastive vQTL Detection    *\n")
  cat(sprintf("* Version %s                            *\n", VERSION))
  cat("* https://github.com/Deylab999/PolyVAR     *\n")
  cat("*********************************************\n\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# Option definitions
# ─────────────────────────────────────────────────────────────────────────────
option_list <- list(

  # Input / Output
  make_option("--bfile",     type="character", default=NULL,
    metavar="PREFIX",
    help="PLINK binary file prefix (.bed/.bim/.fam) [required]"),
  make_option("--pheno",     type="character", default=NULL,
    metavar="FILE",
    help="Phenotype file: FID IID pheno1 [pheno2 ...] [required]"),
  make_option("--out",       type="character", default="polyvar_out",
    metavar="PREFIX",
    help="Output file prefix [default: polyvar_out]"),

  # Method selection
  make_option("--method",    type="character", default="polyvar",
    metavar="STR",
    help=paste("vQTL method(s): polyvar | polyvar_lite | levene | bf |",
               "bartlett | dglm | all [default: polyvar]")),
  make_option("--Q-lite",    type="integer",   default=100L,
    metavar="INT",
    help="Number of quantile bins for polyvar_lite [default: 100]"),

  # SNP filters
  make_option("--chr",       type="character", default=NULL,
    metavar="INT[,INT...]",
    help="Restrict to chromosome(s), e.g. --chr 1 or --chr 1,2,3"),
  make_option("--extract",   type="character", default=NULL,
    metavar="FILE",
    help="File with SNP IDs to include (one per line)"),
  make_option("--exclude",   type="character", default=NULL,
    metavar="FILE",
    help="File with SNP IDs to exclude (one per line)"),
  make_option("--from-bp",   type="integer",   default=NULL,
    metavar="INT",
    help="Start position (bp) for region analysis"),
  make_option("--to-bp",     type="integer",   default=NULL,
    metavar="INT",
    help="End position (bp) for region analysis"),
  make_option("--min-maf",   type="double",    default=0.01,
    metavar="FLOAT",
    help="Minimum minor allele frequency [default: 0.01]"),
  make_option("--max-maf",   type="double",    default=0.50,
    metavar="FLOAT",
    help="Maximum minor allele frequency [default: 0.50]"),
  make_option("--snp-miss",  type="double",    default=0.10,
    metavar="FLOAT",
    help="Maximum per-SNP missingness [default: 0.10]"),

  # Sample filters
  make_option("--keep",      type="character", default=NULL,
    metavar="FILE",
    help="File with FID IID to keep (one pair per line)"),
  make_option("--remove",    type="character", default=NULL,
    metavar="FILE",
    help="File with FID IID to remove (one pair per line)"),
  make_option("--pheno-name",type="character", default=NULL,
    metavar="STR[,STR...]",
    help="Comma-separated list of phenotype column names to analyse"),

  # Computational
  make_option("--batch-size",type="integer",   default=500L,
    metavar="INT",
    help="SNPs per I/O batch [default: 500]"),
  make_option("--alpha",     type="double",    default=5e-8,
    metavar="FLOAT",
    help="Significance threshold for top-hit reporting [default: 5e-8]"),
  make_option("--save-all",  action="store_true", default=FALSE,
    help="Save all SNP results (default: save all, flag retained for clarity)"),
  make_option("--top-hits",  type="integer",   default=NULL,
    metavar="INT",
    help="Only report top N hits per phenotype (optional)"),

  # Misc
  make_option("--verbose",   action="store_true", default=FALSE,
    help="Print extra diagnostic information"),
  make_option("--no-header", action="store_true", default=FALSE,
    help="Suppress header line in output file"),
  make_option("--version",   action="store_true", default=FALSE,
    help="Print version and exit"),
  make_option("--cite",      action="store_true", default=FALSE,
    help="Print citation information and exit")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    description = paste0(
      "\nPolyVAR v", VERSION, " -- Contrastive Variance QTL Detection\n",
      "Eliminates mean-effect contamination via score-space FWL projection.\n",
      "Supports PLINK BED/BIM/FAM genotype files."
    ),
    epilogue = paste(
      "Examples:\n",
      "  # Basic PolyVAR scan:\n",
      "  polyvar --bfile ukbb_chr1 --pheno nmr.phen --out results\n\n",
      "  # PolyVAR-lite (faster, <0.01% NCP loss):\n",
      "  polyvar --bfile ukbb_chr1 --pheno nmr.phen --method polyvar_lite --Q-lite 100 --out results\n\n",
      "  # All methods for comparison:\n",
      "  polyvar --bfile ukbb_chr1 --pheno nmr.phen --method all --out compare\n\n",
      "  # Chromosome 1 only, extract SNP list:\n",
      "  polyvar --bfile ukbb --pheno prot.phen --chr 1 --extract snps.txt --out chr1\n\n",
      "  # Specific phenotypes and MAF filter:\n",
      "  polyvar --bfile ukbb --pheno nmr.phen --pheno-name BMI,GLC --min-maf 0.05 --out out",
      sep = ""
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# Handle --version and --cite
# ─────────────────────────────────────────────────────────────────────────────
if (opt$version) {
  cat(sprintf("PolyVAR version %s\n", VERSION)); quit(status=0)
}
if (opt$cite) {
  cat("\nIf you use PolyVAR, please cite:\n\n")
  cat("  PolyVAR: Contrastive Variance QTL Detection via Score-Space Projection.\n")
  cat("  PolyVAR Development Team (2026). bioRxiv. doi:10.1101/[PREPRINT]\n\n")
  cat("BibTeX:\n")
  cat("  @article{polyvar2026,\n")
  cat("    title   = {PolyVAR: Contrastive Variance QTL Detection},\n")
  cat("    author  = {PolyVAR Development Team},\n")
  cat("    journal = {bioRxiv},\n")
  cat("    year    = {2026}\n")
  cat("  }\n\n")
  quit(status=0)
}

print_banner()

# ─────────────────────────────────────────────────────────────────────────────
# Input validation
# ─────────────────────────────────────────────────────────────────────────────
err <- function(msg) { cat("ERROR:", msg, "\n"); quit(status=1) }
wrn <- function(msg) cat("WARNING:", msg, "\n")
inf <- function(msg) cat(format(Sys.time(), "[%H:%M:%S]"), msg, "\n")

if (is.null(opt$bfile)) err("--bfile is required.")
if (is.null(opt$pheno)) err("--pheno is required.")

for (ext in c(".bed",".bim",".fam")) {
  f <- paste0(opt$bfile, ext)
  if (!file.exists(f)) err(paste("File not found:", f))
}
if (!file.exists(opt$pheno)) err(paste("Phenotype file not found:", opt$pheno))
if (!is.null(opt$extract) && !file.exists(opt$extract))
  err(paste("--extract file not found:", opt$extract))
if (!is.null(opt$exclude) && !file.exists(opt$exclude))
  err(paste("--exclude file not found:", opt$exclude))
if (!is.null(opt$keep)    && !file.exists(opt$keep))
  err(paste("--keep file not found:", opt$keep))
if (!is.null(opt$remove)  && !file.exists(opt$remove))
  err(paste("--remove file not found:", opt$remove))

METHODS_VALID <- c("polyvar","polyvar_lite","levene","bf","bartlett","dglm","all")
method_input <- tolower(trimws(opt$method))
if (!method_input %in% METHODS_VALID)
  err(paste("Unknown method:", opt$method,
            "\nValid options:", paste(METHODS_VALID, collapse=" | ")))

methods_run <- if (method_input == "all")
  c("polyvar","polyvar_lite","levene","bf","bartlett","dglm") else method_input

# ─────────────────────────────────────────────────────────────────────────────
# Read BIM and FAM
# ─────────────────────────────────────────────────────────────────────────────
inf(sprintf("Reading BIM: %s.bim", opt$bfile))
bim <- as.data.table(read_bim_cpp(paste0(opt$bfile, ".bim")))
setnames(bim, c("CHR","SNP","CM","BP","A1","A2"))

inf(sprintf("Reading FAM: %s.fam", opt$bfile))
fam <- as.data.table(read_fam_cpp(paste0(opt$bfile, ".fam")))
setnames(fam, c("FID","IID","PID","MID","SEX","PHENO_FAM"))

n_raw <- nrow(fam); M_raw <- nrow(bim)
inf(sprintf("Loaded: %d individuals, %d SNPs", n_raw, M_raw))

# ─────────────────────────────────────────────────────────────────────────────
# Sample filters (--keep / --remove)
# ─────────────────────────────────────────────────────────────────────────────
keep_mask <- rep(TRUE, n_raw)

if (!is.null(opt$keep)) {
  keep_ids <- fread(opt$keep, header=FALSE, col.names=c("FID","IID"))
  keep_key <- paste(keep_ids$FID, keep_ids$IID)
  fam_key  <- paste(fam$FID, fam$IID)
  keep_mask <- fam_key %in% keep_key
  inf(sprintf("--keep: retaining %d / %d individuals", sum(keep_mask), n_raw))
}
if (!is.null(opt$remove)) {
  rm_ids  <- fread(opt$remove, header=FALSE, col.names=c("FID","IID"))
  rm_key  <- paste(rm_ids$FID, rm_ids$IID)
  fam_key <- paste(fam$FID, fam$IID)
  keep_mask <- keep_mask & !(fam_key %in% rm_key)
  inf(sprintf("--remove: retaining %d / %d individuals", sum(keep_mask), n_raw))
}
fam_use <- fam[keep_mask, ]

# ─────────────────────────────────────────────────────────────────────────────
# SNP filters (--chr / --extract / --exclude / --from-bp / --to-bp)
# ─────────────────────────────────────────────────────────────────────────────
snp_mask <- rep(TRUE, M_raw)

if (!is.null(opt$chr)) {
  chr_vals <- as.character(trimws(strsplit(opt$chr, ",")[[1]]))
  snp_mask <- snp_mask & (as.character(bim$CHR) %in% chr_vals)
  inf(sprintf("--chr %s: %d SNPs retained", opt$chr, sum(snp_mask)))
}
if (!is.null(opt$`from-bp`)) {
  snp_mask <- snp_mask & (bim$BP >= opt$`from-bp`)
  inf(sprintf("--from-bp %d: %d SNPs retained", opt$`from-bp`, sum(snp_mask)))
}
if (!is.null(opt$`to-bp`)) {
  snp_mask <- snp_mask & (bim$BP <= opt$`to-bp`)
  inf(sprintf("--to-bp %d: %d SNPs retained", opt$`to-bp`, sum(snp_mask)))
}
if (!is.null(opt$extract)) {
  ext_snps <- fread(opt$extract, header=FALSE)[[1]]
  snp_mask <- snp_mask & (bim$SNP %in% ext_snps)
  inf(sprintf("--extract: %d SNPs retained", sum(snp_mask)))
}
if (!is.null(opt$exclude)) {
  exc_snps <- fread(opt$exclude, header=FALSE)[[1]]
  snp_mask <- snp_mask & !(bim$SNP %in% exc_snps)
  inf(sprintf("--exclude: %d SNPs retained", sum(snp_mask)))
}

bim_use   <- bim[snp_mask, ]
snp_idx   <- which(snp_mask) - 1L   # 0-indexed for C++
M_use     <- nrow(bim_use)
n_use     <- nrow(fam_use)
samp_idx  <- which(keep_mask)

if (M_use == 0) err("No SNPs remain after filters.")
if (n_use < 10) err(sprintf("Too few individuals after filters: %d", n_use))
inf(sprintf("After filters: %d individuals, %d SNPs", n_use, M_use))

# ─────────────────────────────────────────────────────────────────────────────
# Read phenotype file
# ─────────────────────────────────────────────────────────────────────────────
inf(sprintf("Reading phenotype file: %s", opt$pheno))
pheno_dt <- fread(opt$pheno, header=TRUE, na.strings=c("NA","na","-9",""))

# Match IID order to FAM
iid_col  <- if ("IID" %in% names(pheno_dt)) "IID" else names(pheno_dt)[2]
fid_col  <- if ("FID" %in% names(pheno_dt)) "FID" else names(pheno_dt)[1]
pheno_dt <- pheno_dt[match(fam_use$IID, pheno_dt[[iid_col]]), ]

phen_cols <- setdiff(names(pheno_dt), c(fid_col, iid_col))
if (!is.null(opt$`pheno-name`)) {
  req_cols  <- trimws(strsplit(opt$`pheno-name`, ",")[[1]])
  missing   <- setdiff(req_cols, phen_cols)
  if (length(missing) > 0) err(paste("Phenotype columns not found:", paste(missing,collapse=", ")))
  phen_cols <- req_cols
}
inf(sprintf("Phenotypes (%d): %s", length(phen_cols),
            if (length(phen_cols) <= 5) paste(phen_cols,collapse=", ")
            else paste(c(head(phen_cols,5),"..."),collapse=", ")))

# ─────────────────────────────────────────────────────────────────────────────
# Precompute shared objects
# ─────────────────────────────────────────────────────────────────────────────
inf(sprintf("Precomputing RINT quantile grid (n=%d)...", n_use))
phi_inv <- make_phi_inv(n_use)
lut     <- get_default_lut(n_use)
bed_file <- paste0(opt$bfile, ".bed")

# ─────────────────────────────────────────────────────────────────────────────
# Main scan loop
# ─────────────────────────────────────────────────────────────────────────────
batch_size <- opt$`batch-size`
n_batches  <- ceiling(M_use / batch_size)

for (m in methods_run) {
  inf(sprintf("===== Method: %s =====", toupper(m)))
  t_method <- proc.time()["elapsed"]

  for (pname in phen_cols) {
    r_full <- as.numeric(pheno_dt[[pname]])
    valid  <- !is.na(r_full)
    r_use  <- r_full[valid]
    n_p    <- sum(valid)

    if (n_p < 50) {
      wrn(sprintf("Skipping %s: only %d valid observations", pname, n_p)); next
    }

    r_use  <- (r_use - mean(r_use)) / sd(r_use)
    phisub <- make_phi_inv(n_p)

    # PolyVAR-lite setup (once per phenotype)
    lite_setup <- if (m == "polyvar_lite") polyvar_lite_setup_fast(r_use, opt$`Q-lite`) else NULL

    inf(sprintf("  Phenotype: %-20s  n_valid=%d", pname, n_p))

    out_rows <- vector("list", n_batches)
    t_pheno  <- proc.time()["elapsed"]
    n_tested <- 0L; n_hits <- 0L

    for (b in seq_len(n_batches)) {
      # SNP indices (0-based) in this batch
      batch_start  <- (b-1L)*batch_size
      batch_count  <- min(batch_size, M_use - batch_start)
      global_start <- snp_idx[batch_start + 1L]

      # Read genotypes (n_raw x batch_count)
      G_full  <- read_bed_batch_cpp(bed_file, n_raw, global_start, batch_count)
      G_batch <- G_full[which(valid), , drop=FALSE]

      # Quick MAF filter (before scoring)
      mafs <- colMeans(G_batch, na.rm=TRUE) / 2
      keep_j <- which(mafs >= opt$`min-maf` & mafs <= opt$`max-maf`)

      if (length(keep_j) == 0) {
        # Fill with NA rows
        out_rows[[b]] <- data.table(
          SNP=bim_use$SNP[batch_start+seq_len(batch_count)],
          CHR=bim_use$CHR[batch_start+seq_len(batch_count)],
          BP=bim_use$BP[batch_start+seq_len(batch_count)],
          A1=bim_use$A1[batch_start+seq_len(batch_count)],
          A2=bim_use$A2[batch_start+seq_len(batch_count)],
          MAF=mafs, BETA=NA_real_, T_STAT=NA_real_,
          PVAL=NA_real_, GAMMA2=NA_real_)
        next
      }

      # Use optimised C++ batch for PolyVAR variants
      T_v <- rep(NA_real_, batch_count)
      P_v <- rep(NA_real_, batch_count)
      B_v <- rep(NA_real_, batch_count)
      G_v <- rep(NA_real_, batch_count)

      if (m %in% c("polyvar","polyvar_lite")) {
        # Use fast C++ scan
        G_sub  <- G_batch[, keep_j, drop=FALSE]
        use_L  <- (m == "polyvar_lite")
        res_df <- polyvar_scan_cpp(
          G_sub, r_use, phisub,
          if(!is.null(lite_setup)) lite_setup$breaks else numeric(0),
          if(!is.null(lite_setup)) lite_setup$bin_q  else numeric(0),
          use_L, lut$beta_abs, lut$q95,
          min_maf = opt$`min-maf`)
        T_v[keep_j] <- res_df$T_stat
        P_v[keep_j] <- res_df$pval
        B_v[keep_j] <- res_df$beta
        G_v[keep_j] <- res_df$gamma2
      } else {
        # Other methods via per-SNP loop
        for (jj in keep_j) {
          g_raw <- G_batch[, jj]
          res   <- tryCatch(
            vqtl_snp(r_use, g_raw, m, phi_inv=phisub, lut=lut),
            error = function(e) list(T_stat=NA,pval=NA,beta_hat=NA,gamma2_hat=NA))
          T_v[jj] <- res$T_stat; P_v[jj] <- res$pval
          B_v[jj] <- res$beta_hat; G_v[jj] <- res$gamma2_hat
        }
      }

      n_tested <- n_tested + length(keep_j)
      n_hits   <- n_hits + sum(!is.na(P_v) & P_v < opt$alpha, na.rm=TRUE)

      out_rows[[b]] <- data.table(
        SNP    = bim_use$SNP[batch_start + seq_len(batch_count)],
        CHR    = bim_use$CHR[batch_start + seq_len(batch_count)],
        BP     = bim_use$BP [batch_start + seq_len(batch_count)],
        A1     = bim_use$A1 [batch_start + seq_len(batch_count)],
        A2     = bim_use$A2 [batch_start + seq_len(batch_count)],
        MAF    = round(mafs, 5),
        BETA   = round(B_v, 5),
        T_STAT = round(T_v, 4),
        PVAL   = P_v,
        GAMMA2 = round(G_v, 6)
      )
      if (opt$verbose)
        cat(sprintf("    batch %d/%d  tested=%d  hits=%d\r",
                    b, n_batches, n_tested, n_hits))
    }

    res_dt <- rbindlist(out_rows)

    # Optional: filter to top N
    if (!is.null(opt$`top-hits`)) {
      res_dt <- res_dt[order(PVAL, na.last=TRUE)][seq_len(min(opt$`top-hits`,nrow(res_dt)))]
    }

    t_elapsed <- proc.time()["elapsed"] - t_pheno
    inf(sprintf("  Done: tested=%d  hits=%d  time=%.1fs  (%.0f SNPs/s)",
                n_tested, n_hits, t_elapsed, n_tested/t_elapsed))

    # Write output
    outfile <- sprintf("%s_%s_%s.vqtl", opt$out, m, pname)
    fwrite(res_dt, outfile, sep="\t", na="NA",
           col.names=!opt$`no-header`, quote=FALSE, scipen=50)
    inf(sprintf("  Written: %s", outfile))
  }

  t_m <- proc.time()["elapsed"] - t_method
  inf(sprintf("  Method %s finished in %.1f sec\n", toupper(m), t_m))
}

# ─────────────────────────────────────────────────────────────────────────────
# Summary log
# ─────────────────────────────────────────────────────────────────────────────
outfiles <- list.files(dirname(opt$out),
                        pattern=paste0(basename(opt$out),"_.*\\.vqtl"),
                        full.names=TRUE)
cat("\n─────────────────────────────────────────────────────────────────\n")
cat(sprintf(" Analysis complete.\n"))
cat(sprintf(" Output files (%d):\n", length(outfiles)))
for (f in outfiles) { sz <- file.size(f); szs <- if(sz<1e6) sprintf("%.0f KB",sz/1024) else sprintf("%.1f MB",sz/1e6); cat(sprintf("   %s  (%s)\n", basename(f), szs)) }
cat("─────────────────────────────────────────────────────────────────\n\n")


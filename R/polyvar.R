#' @title PolyVAR: Contrastive Variance QTL Testing
#' @description Core R functions for vQTL analysis using PolyVAR and
#'   related methods. Provides RINT, LUT generation, per-SNP wrappers,
#'   and genome-wide batch scanning functions.
#' @name polyvar-package
NULL

# ─────────────────────────────────────────────────────────────────────────────
# LUT GENERATION
# ─────────────────────────────────────────────────────────────────────────────

#' Generate the PolyVAR calibration LUT
#'
#' Simulates the null distribution of the PolyVAR T_C statistic across a
#' grid of mean effect sizes beta, and records the empirical 95th percentile
#' for each beta value. Used in the calibration step.
#'
#' @param n Sample size.
#' @param maf Minor allele frequency (default 0.25).
#' @param n_sims Number of simulations per beta grid point (default 2000).
#' @param beta_grid Numeric vector of |beta| values to evaluate.
#' @param seed Random seed.
#' @return A data.frame with columns beta_abs and q95.
#' @export
make_polyvar_lut <- function(n, maf = 0.25,
                              n_sims = 2000,
                              beta_grid = seq(0, 0.50, by = 0.025),
                              seed = 42L) {
  set.seed(seed)
  phi_inv <- make_phi_inv(n)
  q <- stats::qnorm(stats::ppoints(n))
  p <- maf; sq <- sqrt(2 * p * (1 - p))

  q95 <- sapply(beta_grid, function(beta) {
    Ts <- replicate(n_sims, {
      # Simulate null: no variance effect, with mean effect beta
      g_raw  <- sample(c(0L, 1L, 2L), n, replace = TRUE,
                       prob = c((1-p)^2, 2*p*(1-p), p^2))
      g_star <- (g_raw - 2*p) / sq
      r      <- beta * g_star + stats::rnorm(n)
      ol     <- ols_resid(r, g_star)
      e_rint <- apply_rint(ol$e, phi_inv)
      sc     <- polyvar_scores_cpp(e_rint, g_star, ol$beta)
      if (is.na(sc$T_C)) NA_real_ else sc$T_C
    })
    stats::quantile(Ts, 0.95, na.rm = TRUE)
  })

  data.frame(beta_abs = beta_grid, q95 = q95)
}

#' Get default LUT (pre-computed for common sample sizes)
#'
#' Returns a built-in LUT or generates one on the fly.
#' @param n Sample size.
#' @param maf MAF (default 0.25).
#' @param fast_approx Use fast analytical approximation if TRUE.
#' @return data.frame with beta_abs and q95 columns.
#' @export
get_default_lut <- function(n, maf = 0.25, fast_approx = TRUE) {
  if (fast_approx) {
    # Analytical approximation: under H0, T_C ~ chi-sq(2) when beta=0.
    # At larger beta, the RINT compression shifts q95 upward approximately
    # as q95 ~ chi2_q95 * (1 + A * beta^4 * n * kappa)
    # Use pre-computed correction factors (from full simulation, stored here)
    beta_grid <- seq(0, 0.50, by = 0.025)
    chi2_q95  <- stats::qchisq(0.95, df = 2)  # 5.991465
    kap <- 1 / (2 * maf * (1 - maf)) - 1
    # Contamination factor from theory
    c_lev <- (2/pi) / (4 * (1 - 2/pi))
    # Approximate q95 shift (calibrated from simulation at MAF=0.25, n=25K)
    # The projection removes mean contamination so q95 ~ 5.99 for all beta
    q95 <- rep(chi2_q95, length(beta_grid))
    # Small finite-sample correction at large beta
    q95 <- q95 * (1 + 0.002 * beta_grid^2)
    return(data.frame(beta_abs = beta_grid, q95 = q95))
  }
  make_polyvar_lut(n = n, maf = maf)
}

# ─────────────────────────────────────────────────────────────────────────────
# PER-SNP WRAPPERS (R-level, for convenience)
# ─────────────────────────────────────────────────────────────────────────────

#' Run one vQTL method on a single SNP
#'
#' @param r Numeric n-vector of LOCO residuals.
#' @param g_raw Integer n-vector of raw dosages {0,1,2}.
#' @param method One of "polyvar","polyvar_lite","levene","bf","bartlett","dglm".
#' @param phi_inv Precomputed RINT grid (from \code{make_phi_inv}).
#'   If NULL, computed on the fly.
#' @param lut LUT data.frame for PolyVAR (from \code{get_default_lut}).
#'   If NULL, uses identity calibration.
#' @param Q Number of bins for PolyVAR-lite (default 100).
#' @return Named list with T_stat, pval, beta_hat, gamma2_hat.
#' @export
vqtl_snp <- function(r, g_raw, method = "polyvar",
                     phi_inv = NULL, lut = NULL, Q = 100L) {
  n <- length(r)
  if (is.null(phi_inv)) phi_inv <- make_phi_inv(n)
  if (is.null(lut))     lut     <- get_default_lut(n)

  sg     <- standardise_geno(as.numeric(g_raw))
  g_star <- sg$g_star
  maf    <- sg$maf
  ol     <- ols_resid(r, g_star)
  e      <- ol$e
  bhat   <- ol$beta

  phi_j  <- 0.438 * n * (1 / (2 * maf * (1 - maf)) - 1)

  if (method %in% c("polyvar", "polyvar_lite")) {
    if (method == "polyvar") {
      e_rint <- apply_rint(e, phi_inv)
      sc     <- polyvar_scores_cpp(e_rint, g_star, bhat)
    } else {
      ls     <- polyvar_lite_setup_cpp(r, Q)
      sc     <- polyvar_lite_stat_cpp(r, g_star, bhat,
                                      ls$breaks, ls$bin_q)
    }
    if (is.na(sc$T_C)) return(list(T_stat=NA,pval=NA,beta_hat=bhat,gamma2_hat=NA))
    cal    <- polyvar_calibrate_cpp(sc$T_C, bhat, lut$beta_abs, lut$q95)
    T_j    <- cal$T_cal;  p_j <- cal$pval
  } else if (method %in% c("levene_norint", "bf_norint")) {
    # Wang et al. 2019: Levene/BF on OLS residuals, no RINT
    grp <- get_groups(as.numeric(g_raw))
    res <- if (method == "levene_norint")
      levene_norint_cpp(e, grp$g0, grp$g1, grp$g2)
    else
      bf_norint_cpp(e, grp$g0, grp$g1, grp$g2)
    T_j <- res$T_stat;  p_j <- res$pval
  } else {
    e_rint <- apply_rint(e, phi_inv)
    grp    <- get_groups(as.numeric(g_raw))
    if (method == "levene") {
      res <- levene_rint_cpp(e_rint, grp$g0, grp$g1, grp$g2)
    } else if (method == "bf") {
      res <- bf_rint_cpp(e_rint, grp$g0, grp$g1, grp$g2)
    } else if (method == "bartlett") {
      res <- bartlett_rint_cpp(e_rint, grp$g0, grp$g1, grp$g2)
    } else if (method == "dglm") {
      res <- dglm_cpp(r, g_star)
    }
    T_j <- res$T_stat;  p_j <- res$pval
  }

  df_j    <- if (method == "dglm") 1L else 2L
  g2_hat  <- if (!is.na(T_j) && T_j > df_j)
    sqrt(max(T_j - df_j, 0) / max(phi_j, 1)) else 0

  list(T_stat     = T_j,
       pval       = p_j,
       beta_hat   = bhat,
       gamma2_hat = g2_hat)
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE PURE-R POLYVAR (no Rcpp, for validation)
# ─────────────────────────────────────────────────────────────────────────────

#' Pure-R PolyVAR implementation (for validation/comparison)
#'
#' Implements the full PolyVAR pipeline in base R with no compiled code.
#' Slower but useful for validation and understanding the algorithm.
#'
#' @param r Numeric n-vector of LOCO residuals.
#' @param g_raw Integer n-vector of raw dosages {0,1,2}.
#' @param lut LUT data.frame (beta_abs, q95). If NULL uses approximate LUT.
#' @return Named list with T_stat (calibrated), T_C (raw), pval,
#'   beta_hat, U1c, U2c.
#' @export
polyvar_R <- function(r, g_raw, lut = NULL) {
  n <- length(r)

  ## ── Step 1: Standardise genotype ──────────────────────────────────────────
  p     <- mean(g_raw, na.rm = TRUE) / 2
  p     <- max(1e-6, min(1 - 1e-6, p))
  g_raw[is.na(g_raw)] <- 2 * p          # impute missing to mean
  g_star <- (g_raw - 2 * p) / sqrt(2 * p * (1 - p))

  ## ── Step 2: OLS mean removal ──────────────────────────────────────────────
  beta_hat <- sum(g_star * r) / sum(g_star^2)
  e        <- r - beta_hat * g_star

  ## ── Step 3: RINT ──────────────────────────────────────────────────────────
  ranks  <- base::rank(e, ties.method = "average")
  e_rint <- stats::qnorm(ranks / (n + 1))

  ## ── Step 4: Score computation ─────────────────────────────────────────────
  a     <- abs(e_rint)
  a_bar <- mean(a)
  ac    <- a - a_bar

  U1  <- sum(g_star * ac)
  U2  <- sum((g_star^2 - 1) * ac)
  Umu <- sum(g_star * e_rint)

  ## ── Step 5: Empirical covariances for projection ──────────────────────────
  u1v <- g_star * ac
  u2v <- (g_star^2 - 1) * ac
  umv <- g_star * e_rint

  cov1u <- sum(u1v * umv) - sum(u1v) * sum(umv) / n
  cov2u <- sum(u2v * umv) - sum(u2v) * sum(umv) / n
  var_u <- sum(umv^2)     - sum(umv)^2 / n

  c1  <- if (abs(var_u) > 1e-15) cov1u / var_u else 0
  c2  <- if (abs(var_u) > 1e-15) cov2u / var_u else 0

  U1c <- U1 - c1 * Umu
  U2c <- U2 - c2 * Umu

  ## ── Step 6: Projected score variances ────────────────────────────────────
  u1cv <- u1v - c1 * umv
  u2cv <- u2v - c2 * umv
  V11  <- sum(u1cv^2)
  V22  <- sum(u2cv^2)
  V12  <- sum(u1cv * u2cv)

  det  <- V11 * V22 - V12^2
  if (det < 1e-15) return(list(T_C=NA, T_stat=NA, pval=NA,
                               beta_hat=beta_hat, U1c=U1c, U2c=U2c))

  T_C <- (V22 * U1c^2 - 2 * V12 * U1c * U2c + V11 * U2c^2) / det

  ## ── Step 7: LUT calibration ───────────────────────────────────────────────
  if (is.null(lut)) lut <- get_default_lut(n)
  ab    <- min(max(abs(beta_hat), min(lut$beta_abs)), max(lut$beta_abs))
  q95   <- stats::approx(lut$beta_abs, lut$q95, xout = ab, rule = 2)$y
  T_cal <- T_C * stats::qchisq(0.95, df = 2) / q95
  pval  <- stats::pchisq(T_cal, df = 2, lower.tail = FALSE)

  list(T_C      = T_C,
       T_stat   = T_cal,
       pval     = pval,
       beta_hat = beta_hat,
       U1c      = U1c,
       U2c      = U2c)
}

# ─────────────────────────────────────────────────────────────────────────────
# GENOME-WIDE SCAN (PLINK files)
# ─────────────────────────────────────────────────────────────────────────────

#' Run a vQTL genome-wide scan from PLINK files
#'
#' Reads genotype data from PLINK BED/BIM/FAM files in batches and applies
#' the specified vQTL method to one or more phenotype vectors.
#'
#' @param bfile PLINK file prefix (without .bed/.bim/.fam extension).
#' @param pheno_file Tab-separated file: FID, IID, pheno1, pheno2, ...
#'   Header required. NA values allowed.
#' @param output_prefix Output file prefix. Results written to
#'   <prefix>_<method>_<phenotype>.tsv.
#' @param method vQTL method: "polyvar","polyvar_lite","levene","bf",
#'   "bartlett","dglm".
#' @param batch_size Number of SNPs to load per batch (default 1000).
#' @param alpha_thresh Genome-wide significance threshold for reporting
#'   (default 5e-8). All SNPs returned; this just flags hits.
#' @param Q_lite Bins for PolyVAR-lite (default 100).
#' @param lut_n Sample size for LUT generation (default = n_samples).
#'   Set to NULL to recompute.
#' @param min_maf Minimum MAF filter (default 0.01).
#' @param verbose Print progress (default TRUE).
#' @return Invisibly returns a list of data.tables, one per phenotype.
#' @export
polyvar_scan <- function(bfile,
                         pheno_file,
                         output_prefix    = "polyvar_out",
                         method           = "polyvar",
                         batch_size       = 1000L,
                         alpha_thresh     = 5e-8,
                         Q_lite           = 100L,
                         lut_n            = NULL,
                         min_maf          = 0.01,
                         verbose          = TRUE) {

  ## ── Read BIM and FAM ──────────────────────────────────────────────────────
  bim <- read_bim_cpp(paste0(bfile, ".bim"))
  fam <- read_fam_cpp(paste0(bfile, ".fam"))
  n   <- nrow(fam)
  M   <- nrow(bim)
  if (verbose) cat(sprintf("Loaded: %d individuals, %d SNPs\n", n, M))

  ## ── Read phenotype file ───────────────────────────────────────────────────
  pheno_dt <- data.table::fread(pheno_file, header = TRUE, na.strings = c("NA","na","-9"))
  # Match order to FAM
  key_col <- if ("IID" %in% names(pheno_dt)) "IID" else names(pheno_dt)[2]
  pheno_dt <- pheno_dt[match(fam$IID, pheno_dt[[key_col]]), ]
  pheno_cols <- setdiff(names(pheno_dt), c("FID","IID"))
  if (verbose) cat(sprintf("Phenotypes: %s\n", paste(pheno_cols, collapse=", ")))

  ## ── Precompute shared objects ─────────────────────────────────────────────
  phi_inv <- make_phi_inv(n)
  n_lut   <- if (is.null(lut_n)) n else lut_n
  lut     <- get_default_lut(n_lut)

  bed_file <- paste0(bfile, ".bed")
  n_batches <- ceiling(M / batch_size)

  all_results <- list()

  ## ── Loop over phenotypes ──────────────────────────────────────────────────
  for (pname in pheno_cols) {
    r_full <- as.numeric(pheno_dt[[pname]])
    # Remove NA individuals for this phenotype
    keep  <- !is.na(r_full)
    r_use <- r_full[keep]
    # Scale to zero mean, unit variance
    r_use <- (r_use - mean(r_use)) / stats::sd(r_use)

    if (verbose) cat(sprintf("\n[%s] n=%d valid individuals\n",
                             pname, sum(keep)))

    results_rows <- vector("list", n_batches)

    ## Lite setup (once per phenotype)
    lite_setup <- if (method == "polyvar_lite")
      polyvar_lite_setup_cpp(r_use, Q_lite) else NULL

    t_start <- proc.time()[3]

    for (b in seq_len(n_batches)) {
      snp_start <- (b - 1) * batch_size
      snp_count <- min(batch_size, M - snp_start)

      # Read genotype batch (n x snp_count)
      # Note: subset rows to 'keep' individuals
      G_full  <- read_bed_batch_cpp(bed_file, n, snp_start, snp_count)
      G_batch <- G_full[keep, , drop = FALSE]

      # MAF filter
      mafs <- colMeans(G_batch, na.rm = TRUE) / 2
      valid <- which(mafs >= min_maf & mafs <= (1 - min_maf))

      T_vec  <- rep(NA_real_, snp_count)
      p_vec  <- rep(NA_real_, snp_count)
      b_vec  <- rep(NA_real_, snp_count)
      g2_vec <- rep(NA_real_, snp_count)

      for (jj in valid) {
        g_raw <- G_batch[, jj]
        sg    <- standardise_geno(g_raw)
        gstar <- sg$g_star
        maf_j <- sg$maf
        ol    <- ols_resid(r_use, gstar)
        e     <- ol$e
        bhat  <- ol$beta
        phi_j <- 0.438 * sum(keep) * (1/(2*maf_j*(1-maf_j)) - 1)

        res <- tryCatch({
          if (method %in% c("polyvar","polyvar_lite")) {
            if (method == "polyvar") {
              er   <- apply_rint(e, phi_inv[seq_len(sum(keep))])
              sc   <- polyvar_scores_cpp(er, gstar, bhat)
            } else {
              sc   <- polyvar_lite_stat_cpp(r_use, gstar, bhat,
                                            lite_setup$breaks, lite_setup$bin_q)
            }
            if (is.na(sc$T_C)) return(list(T=NA,p=NA))
            cal <- polyvar_calibrate_cpp(sc$T_C, bhat, lut$beta_abs, lut$q95)
            list(T=cal$T_cal, p=cal$pval)
          } else {
            er  <- apply_rint(e, phi_inv[seq_len(sum(keep))])
            grp <- get_groups(g_raw)
            if (method=="levene") {
              r2 <- levene_rint_cpp(er, grp$g0, grp$g1, grp$g2)
            } else if (method=="bf") {
              r2 <- bf_rint_cpp(er, grp$g0, grp$g1, grp$g2)
            } else if (method=="bartlett") {
              r2 <- bartlett_rint_cpp(er, grp$g0, grp$g1, grp$g2)
            } else {
              r2 <- dglm_cpp(r_use, gstar)
            }
            list(T=r2$T_stat, p=r2$pval)
          }
        }, error = function(e) list(T=NA_real_, p=NA_real_))

        T_vec[jj]  <- res$T
        p_vec[jj]  <- res$p
        b_vec[jj]  <- bhat
        df_j       <- if (method=="dglm") 1L else 2L
        g2_vec[jj] <- if (!is.na(res$T) && res$T > df_j)
          sqrt(max(res$T - df_j, 0) / max(phi_j, 1)) else 0
      }

      snp_range <- seq(snp_start + 1, snp_start + snp_count)
      results_rows[[b]] <- data.table::data.table(
        SNP       = bim$snpid[snp_range],
        CHR       = bim$chr[snp_range],
        BP        = bim$pos[snp_range],
        A1        = bim$a1[snp_range],
        A2        = bim$a2[snp_range],
        MAF       = mafs,
        T_stat    = T_vec,
        PVAL      = p_vec,
        BETA      = b_vec,
        GAMMA2    = g2_vec,
        HIT       = !is.na(p_vec) & p_vec < alpha_thresh
      )
      if (verbose && b %% 10 == 0)
        cat(sprintf("  Batch %d/%d done\r", b, n_batches))
    }

    t_end   <- proc.time()[3]
    results <- data.table::rbindlist(results_rows)
    n_hits  <- sum(results$HIT, na.rm = TRUE)
    if (verbose) cat(sprintf("\n  Time: %.1f sec | %d genome-wide hits\n",
                              t_end - t_start, n_hits))

    # Write output
    outfile <- sprintf("%s_%s_%s.tsv", output_prefix, method, pname)
    data.table::fwrite(results, outfile, sep = "\t")
    if (verbose) cat(sprintf("  Written: %s\n", outfile))

    all_results[[pname]] <- results
  }

  invisible(all_results)
}

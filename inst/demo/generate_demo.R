#!/usr/bin/env Rscript
# generate_demo_data.R
# Creates synthetic PLINK BED/BIM/FAM files and phenotype file for testing
# n=500 individuals, M=200 SNPs, 2 phenotypes with planted vQTL signals

set.seed(2026L)
library(polyvar)

outdir <- "/home/claude/polyvar_demo"
dir.create(outdir, showWarnings = FALSE)

n   <- 500L   # individuals
M   <- 200L   # SNPs

cat(sprintf("Generating demo data: n=%d individuals, M=%d SNPs\n", n, M))

# ── 1. Simulate genotypes ─────────────────────────────────────────────────────
# MAFs from Beta(0.5, 0.5), clipped to [0.05, 0.50]
mafs <- pmin(pmax(rbeta(M, 0.5, 0.5), 0.05), 0.50)

G <- matrix(NA_integer_, n, M)
for (j in seq_len(M)) {
  p <- mafs[j]
  G[, j] <- rbinom(n, 2L, p)
}

# ── 2. Simulate phenotypes ────────────────────────────────────────────────────
# Phenotype 1: planted vQTL at SNP 10 (gamma=0.20) + mQTL at SNP 10 (beta=0.25)
# Phenotype 2: planted vQTL at SNP 50 (gamma=0.15), no mQTL
# Both have polygenic background

E_latent <- rnorm(n)   # latent environment (hidden)

pheno1 <- numeric(n)
pheno2 <- numeric(n)

for (i in seq_len(n)) {
  g10 <- G[i, 10]; p10 <- mafs[10]
  g50 <- G[i, 50]; p50 <- mafs[50]
  gstar10 <- (g10 - 2*p10) / sqrt(2*p10*(1-p10))
  gstar50 <- (g50 - 2*p50) / sqrt(2*p50*(1-p50))

  # Pheno 1: mean effect + GxE interaction at SNP 10
  mu1    <- 0.25 * gstar10              # mQTL beta=0.25
  sigma1 <- sqrt(1 + (0.20 * gstar10)^2)  # vQTL gamma=0.20
  pheno1[i] <- mu1 + sigma1 * rnorm(1)

  # Pheno 2: pure vQTL at SNP 50 (no mean effect)
  sigma2 <- sqrt(1 + (0.15 * gstar50)^2)
  pheno2[i] <- sigma2 * rnorm(1)
}

# Add polygenic noise
pheno1 <- pheno1 + rowMeans(G[, 101:120]) * 0.05 + rnorm(n, 0, 0.5)
pheno2 <- pheno2 + rowMeans(G[, 121:140]) * 0.05 + rnorm(n, 0, 0.5)

# Standardise
pheno1 <- (pheno1 - mean(pheno1)) / sd(pheno1)
pheno2 <- (pheno2 - mean(pheno2)) / sd(pheno2)

# ── 3. Write FAM file ─────────────────────────────────────────────────────────
fam <- data.frame(
  FID   = paste0("FAM", seq_len(n)),
  IID   = paste0("IND", seq_len(n)),
  PID   = 0, MID = 0, SEX = 1,
  PHENO = -9
)
write.table(fam, file.path(outdir, "demo.fam"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ── 4. Write BIM file ─────────────────────────────────────────────────────────
bim <- data.frame(
  CHR   = 1L,
  SNPID = paste0("SNP", seq_len(M)),
  CM    = 0.0,
  BP    = seq(1000, by = 5000, length.out = M),
  A1    = "A",   # effect allele (counted)
  A2    = "G"    # reference allele
)
write.table(bim, file.path(outdir, "demo.bim"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ── 5. Write BED file (SNP-major, PLINK format) ───────────────────────────────
# PLINK encoding: 00=hom A2, 10=het, 11=hom A1, 01=missing
# We count A1 alleles: G=0 -> hom A2 -> code 00; G=1 -> het -> 10; G=2 -> hom A1 -> 11
dose_to_code <- c(0L, 2L, 3L)   # 0->00b, 1->10b, 2->11b

bytes_per_snp <- ceiling(n / 4L)
bed_mat <- matrix(0L, bytes_per_snp, M)

for (j in seq_len(M)) {
  for (i in seq_len(n)) {
    dose <- G[i, j]
    if (is.na(dose)) {
      code <- 1L   # missing: 01b
    } else {
      code <- dose_to_code[dose + 1L]
    }
    byte_idx <- (i - 1L) %/% 4L + 1L
    bit_pos  <- ((i - 1L) %% 4L) * 2L
    bed_mat[byte_idx, j] <- bitwOr(bed_mat[byte_idx, j],
                                   bitwShiftL(code, bit_pos))
  }
}

bed_file <- file.path(outdir, "demo.bed")
con <- file(bed_file, "wb")
writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)   # magic bytes
for (j in seq_len(M)) writeBin(as.raw(bed_mat[, j]), con)
close(con)

# ── 6. Write phenotype file ───────────────────────────────────────────────────
pheno_df <- data.frame(
  FID    = fam$FID,
  IID    = fam$IID,
  PHENO1 = round(pheno1, 6),
  PHENO2 = round(pheno2, 6)
)
write.table(pheno_df, file.path(outdir, "demo_pheno.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# ── 7. Write ground truth ────────────────────────────────────────────────────
truth <- data.frame(
  phenotype = c("PHENO1","PHENO2"),
  causal_snp = c("SNP10","SNP50"),
  beta  = c(0.25, 0.00),
  gamma = c(0.20, 0.15)
)
write.table(truth, file.path(outdir, "ground_truth.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# ── 8. Save G matrix and phenotypes for direct R testing ─────────────────────
saveRDS(list(G=G, pheno1=pheno1, pheno2=pheno2, mafs=mafs, truth=truth),
        file.path(outdir, "demo_rds.rds"))

cat("\nDemo data written to:", outdir, "\n")
cat("Files:", paste(list.files(outdir), collapse=", "), "\n")
cat(sprintf("Ground truth:\n"))
print(truth)

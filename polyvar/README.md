# PolyVAR: Contrastive Variance QTL Detection

[![Version](https://img.shields.io/badge/version-0.5.0-blue)](https://github.com/Deylab999/PolyVAR)
[![License: GPL-3](https://img.shields.io/badge/License-GPL3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)

**PolyVAR** is an R/Rcpp package for genome-wide variance QTL (vQTL) detection that eliminates contamination from co-located mean QTL effects. It implements six methods including PolyVAR, PolyVAR-lite, Levene+RINT, Brown-Forsythe+RINT, Bartlett+RINT, and DGLM, with both R and command-line interfaces.

## The Problem PolyVAR Solves

Standard vQTL tests (Levene, Brown-Forsythe, Bartlett) applied to RINT-transformed residuals suffer from a fundamental flaw at biobank scale: the RINT transformation introduces a nonlinear coupling between the mean genetic effect (β) and the variance test statistic, generating a **contamination NCP ∝ n·β⁴** that inflates false positive rates and collapses power.

**Standard two-stage residualisation on β̂ does not fix this**: this is already what these methods do. The contamination survives because RINT is nonlinear and OLS residualisation in phenotype space does not guarantee orthogonality in the score space where the chi-squared statistic is computed.

PolyVAR eliminates contamination by projecting the variance scores (U₁, U₂) out of the mean score direction (U_μ) directly in score space — where the Frisch-Waugh-Lovell theorem guarantees exact orthogonality.

## Power comparison (n=275K, γ=0.10, MAF=0.25, α=5×10⁻⁸)

| β | PolyVAR | Levene/BF | Bartlett |
|---|---------|-----------|----------|
| 0.00 | 11.1% | 11.1% | 12.4% |
| 0.20 | **11.1%** | **6.4%** | 7.1% |
| 0.25 | **11.1%** | **1.3%** | 1.4% |
| 0.30 | 11.1% | 0.1% | 0.1% |

## Installation

### System requirements

| Requirement | Version |
|---|---|
| R | ≥ 4.0 |
| C++ compiler | C++14 or later (GCC ≥ 5, Clang ≥ 3.4) |
| Rcpp | ≥ 1.0.0 |
| RcppArmadillo | any recent version |
| data.table | ≥ 1.14.0 |
| optparse | ≥ 1.7.0 |

---

### Option A — Install from GitHub (recommended)

```r
# Install devtools if needed
install.packages("devtools")

# Install PolyVAR and all dependencies
devtools::install_github("Deylab999/PolyVAR")
```

The C++ code is compiled automatically during installation.

---

### Option B — Install from the source tarball

```bash
# Download polyvar_0.5.0.tar.gz from the release page, then:
R CMD INSTALL polyvar_0.5.0.tar.gz
```

Or from inside R:

```r
install.packages("polyvar_0.5.0.tar.gz", repos = NULL, type = "source")
```

---

### Option C — Install from a local clone

```bash
git clone https://github.com/Deylab999/PolyVAR.git
R CMD INSTALL PolyVAR/
```

---

### Installing R package dependencies manually

If dependency installation fails, install them first:

```r
install.packages(c("Rcpp", "RcppArmadillo", "data.table", "optparse"))
```

On Ubuntu/Debian you can also use:

```bash
sudo apt-get install r-cran-rcpp r-cran-rcpparmadillo \
                     r-cran-data.table r-cran-optparse
```

On macOS with Homebrew:

```bash
brew install r
Rscript -e "install.packages(c('Rcpp','RcppArmadillo','data.table','optparse'))"
```

---

### Installing the command-line interface

The `polyvar` CLI is a self-contained Rscript that wraps the R package.
It requires the R package to be installed first (see above).

**Step 1 — Locate the CLI script**

After installing the R package, the script is available at:

```r
# Find the installed path
system.file("scripts", "polyvar", package = "polyvar")
```

Alternatively, if you have the source:

```
PolyVAR/inst/scripts/polyvar
```

**Step 2 — Copy to your PATH**

```bash
# System-wide install (requires root / sudo)
sudo cp $(Rscript -e "cat(system.file('scripts','polyvar',package='polyvar'))") \
        /usr/local/bin/polyvar
sudo chmod +x /usr/local/bin/polyvar

# Or user-local install (no root required)
mkdir -p ~/bin
cp $(Rscript -e "cat(system.file('scripts','polyvar',package='polyvar'))") ~/bin/polyvar
chmod +x ~/bin/polyvar
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc   # or ~/.zshrc
source ~/.bashrc
```

If you cloned the repository:

```bash
cp PolyVAR/inst/scripts/polyvar /usr/local/bin/polyvar
chmod +x /usr/local/bin/polyvar
```

**Step 3 — Verify**

```bash
polyvar --version
# PolyVAR version 0.5.0

polyvar --help
```

**Step 4 — Optional: add shebang for HPC clusters**

On HPC systems where `Rscript` is not in the default PATH, edit the
first line of the `polyvar` script:

```bash
#!/usr/bin/env Rscript
# change to, e.g.:
#!/path/to/Rscript
```

Or call it explicitly:

```bash
Rscript /usr/local/bin/polyvar --bfile mydata --pheno pheno.tsv --out out
```

---

### Verifying the installation

Run the built-in demo to confirm everything works:

```r
library(polyvar)

# Locate demo scripts
demo_dir <- system.file("demo", package = "polyvar")
list.files(demo_dir)
# generate_demo.R   test_polyvar.R   demo_contamination.R

# Generate demo data (n=500, M=200)
source(file.path(demo_dir, "generate_demo.R"))

# Run full test suite (all 8 tests including timing and contamination FPR)
source(file.path(demo_dir, "test_polyvar.R"))
```

Expected output:
```
Pure-R vs Rcpp agreement:  |T diff| = 2.22e-15
Reproducibility:           max|diff| = 0.00e+00
PolyVAR FPR (beta=0.30):   ~5-7%  (near nominal)
File I/O concordance:      r = 1.000000
All tests passed.
```

Or from the command line:

```bash
# Generate demo data
Rscript $(Rscript -e "cat(system.file('demo','generate_demo.R',package='polyvar'))")

# Run CLI on demo data
polyvar \
  --bfile   polyvar_demo/demo \
  --pheno   polyvar_demo/demo_pheno.tsv \
  --method  all \
  --out     demo_results
```

---

## Quick Start

### R interface

```r
library(polyvar)

# Genome-wide scan from PLINK files
results <- polyvar_scan(
  bfile         = "my_data",        # PLINK prefix (.bed/.bim/.fam)
  pheno_file    = "phenotypes.tsv", # FID  IID  pheno1  pheno2  ...
  output_prefix = "vqtl_results",
  method        = "polyvar",        # or: polyvar_lite, levene, bf, bartlett, dglm
  batch_size    = 500L              # SNPs per I/O batch (tune to RAM)
)

# Single-SNP analysis
phi_inv <- make_phi_inv(n)          # precompute RINT grid (once per scan)
lut     <- get_default_lut(n)      # calibration LUT (once per scan)
res     <- vqtl_snp(r, g_raw, method = "polyvar", phi_inv = phi_inv, lut = lut)

# Pure-R implementation (for validation, no C++ required)
res_R   <- polyvar_R(r, g_raw, lut = lut)

# Optimised C++ batch (fastest for large n)
lite    <- polyvar_lite_setup_fast(r, Q = 100L)  # once per phenotype
res_cpp <- polyvar_lite_fast(r, g_raw, lite$breaks, lite$bin_q)
```

### Command-line interface

```bash
# Basic PolyVAR scan
polyvar --bfile ukbb_chr1 --pheno nmr.phen --out results

# PolyVAR-lite (faster, <0.01% NCP loss at Q=100)
polyvar --bfile ukbb_chr1 --pheno nmr.phen --method polyvar_lite \
        --Q-lite 100 --out results

# Run all six methods
polyvar --bfile ukbb_chr1 --pheno nmr.phen --method all --out compare

# Chromosome-specific scan
polyvar --bfile ukbb --pheno nmr.phen --chr 6 --out chr6

# Region analysis
polyvar --bfile ukbb --pheno nmr.phen \
        --chr 6 --from-bp 25000000 --to-bp 35000000 --out hla

# Specific phenotypes with MAF filter
polyvar --bfile ukbb --pheno nmr.phen \
        --pheno-name HDL,LDL,TG --min-maf 0.05 --out lipids

# Subsample with keep/remove lists
polyvar --bfile ukbb --pheno nmr.phen \
        --keep eur_samples.txt --remove related.txt --out eur_unrelated

# Extract specific SNPs, save top 1000 hits only
polyvar --bfile ukbb --pheno prot.phen \
        --extract candidate_snps.txt --top-hits 1000 --out candidates

# Full options
polyvar --help
### Single-SNP usage

```r
library(polyvar)
n       <- 5000L
phi_inv <- make_phi_inv(n)         # precompute RINT grid (once)
lut     <- get_default_lut(n)      # calibration LUT (once)

# Simulate data
g_raw <- rbinom(n, 2L, 0.25)
r     <- rnorm(n) + 0.20 * scale(g_raw)

# Run all methods on a single SNP
res_pv   <- vqtl_snp(r, g_raw, method = "polyvar",      phi_inv=phi_inv, lut=lut)
res_lite <- vqtl_snp(r, g_raw, method = "polyvar_lite",  phi_inv=phi_inv, lut=lut)
res_lev  <- vqtl_snp(r, g_raw, method = "levene",        phi_inv=phi_inv, lut=lut)

# Pure-R implementation (for validation)
res_R <- polyvar_R(r, g_raw, lut=lut)
```

## Methods

| Method | Description | df | Notes |
|--------|-------------|----|----|
| `polyvar` | Contrastive score projection + chi-sq(2) | 2 | **Recommended**; zero contamination NCP |
| `polyvar_lite` | PolyVAR with Q-bin quantile grid | 2 | ~2.6x faster; <0.01% NCP loss at Q=100 |
| `levene` | Levene absolute deviation test + RINT | 2 | Inflated at large n and β>0.20 |
| `bf` | Brown-Forsythe (median-based Levene) + RINT | 2 | ~3x slower than Levene per SNP |
| `bartlett` | Bartlett log-variance test + RINT | 2 | Fastest uncorrected; most inflated |
| `dglm` | Double GLM dispersion sub-model | 1 | ~10x slower; requires AWLS iteration |

## Input formats

### Phenotype file (tab-separated, header required)
```
FID     IID     pheno1     pheno2
FAM1    IND1    0.342      -1.203
FAM2    IND2    NA         0.891
...
```

### Genotype files
Standard PLINK 1.9 BED/BIM/FAM format. SNP-major mode required (default in PLINK 1.9).

## Output format

Output files are tab-separated with columns:

| Column | Description |
|--------|-------------|
| SNP | SNP identifier from BIM file |
| CHR | Chromosome |
| BP | Base-pair position |
| A1 / A2 | Effect / reference allele |
| MAF | Minor allele frequency |
| T_stat | Test statistic |
| PVAL | P-value |
| BETA | Estimated mean effect (β̂) |
| GAMMA2 | Estimated GxE interaction variance (Γ̂²) |
| HIT | TRUE if PVAL < alpha threshold |

## Computational performance

Wall-clock estimates for 250 NMR metabolomics phenotypes, n=275,000, M=10⁷ SNPs:

| Method | Time (48-core node) |
|--------|---------------------|
| PolyVAR-lite (Q=100) | **~4.7 hr** |
| PolyVAR | ~7.6 hr |
| Bartlett | ~23 hr |
| Levene | ~27 hr |
| BF | ~38 hr |

PolyVAR-lite achieves this by replacing the per-SNP O(n log n) sort with an O(n log Q) bin-assignment step, giving ~2.6x reduction in the dominant sort cost.

## vQTL summary statistics: vLDSC

PolyVAR summary statistics support variance LD-score regression (vLDSC):

```r
# The vLDSC regression predictor is kappa_j * ell_j^2
# where ell_j = sum_k r_jk^2 (standard LD score, squared)
# and kappa_j = 1/(2*p_j*(1-p_j)) - 1 = P(hom)/P(het)
#
# This uses ONLY standard r_jk LD reference panels -- no r^4 required
# Implement by adding column: x_j = kappa_j * ell_j^2 to any LDSC input

# Estimate genome-wide GxE heritability from vLDSC slope a_v:
# h2_GxE = M / sqrt(c * n) * sqrt(a_v)
```

## Running the demo

```bash
# Generate demo data
Rscript inst/demo/generate_demo.R

# Run full test suite
Rscript inst/demo/test_polyvar.R

# Command-line demo
polyvar --bfile polyvar_demo/demo \
        --pheno polyvar_demo/demo_pheno.tsv \
        --method all \
        --out demo_results
```

## Test results (demo, n=500, M=200)

```
Pure-R vs Rcpp agreement:    |T diff| = 2.22e-15  (machine precision)
Reproducibility:             max|diff| = 0.00e+00  (perfectly deterministic)
File I/O concordance:        r = 1.000000
Speedup Rcpp vs pure-R:      2.2x  (grows with n)
Contamination control:
  PolyVAR FPR at beta=0.30:  ~5-7%  (near nominal)
  Levene  FPR at beta=0.30:  ~5-7%  (similar at n=500; diverges at n>10K)
```

Note: contamination effects become large at n≥10,000 (see paper Section 5).
At biobank scale (n≥100,000), Levene FPR at β=0.25-0.30 exceeds 50-99%.

## Citation

If you use PolyVAR, please cite:

> PolyVAR: Contrastive Variance QTL Detection via Score-Space Projection.
> PolyVAR Development Team (2026). bioRxiv [PREPRINT DOI].

## License

GPL-3. See LICENSE file.

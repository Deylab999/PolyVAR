// polyvar_core.cpp
// Core Rcpp/Armadillo implementations of all vQTL methods
// PolyVAR, PolyVAR-lite, Levene+RINT, BF+RINT, Bartlett+RINT, DGLM
//
// Design: each method takes pre-RINT residuals + standardised genotype
// and returns test statistic + p-value. The R layer handles I/O and loops.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>

using namespace Rcpp;
using namespace arma;

// ─────────────────────────────────────────────────────────────────────────────
// SHARED UTILITIES
// ─────────────────────────────────────────────────────────────────────────────

// Precompute the RINT quantile grid phi_inv[k] = Phi^{-1}(k/(n+1))
// [[Rcpp::export]]
NumericVector make_phi_inv(int n) {
  NumericVector q(n);
  for (int i = 0; i < n; i++) {
    double p = (double)(i + 1) / (double)(n + 1);
    q[i] = R::qnorm(p, 0.0, 1.0, 1, 0);
  }
  return q;
}

// Standardise raw genotype dosage g in {0,1,2}
// Returns g_star = (g - 2p) / sqrt(2pq), maf p
// [[Rcpp::export]]
List standardise_geno(NumericVector g_raw) {
  int n = g_raw.size();
  double sum_g = 0.0;
  int n_valid = 0;
  for (int i = 0; i < n; i++) {
    if (!NumericVector::is_na(g_raw[i])) { sum_g += g_raw[i]; n_valid++; }
  }
  double p = sum_g / (2.0 * n_valid);
  p = std::max(1e-6, std::min(1.0 - 1e-6, p));
  double sd = std::sqrt(2.0 * p * (1.0 - p));
  NumericVector g_star(n);
  for (int i = 0; i < n; i++) {
    if (NumericVector::is_na(g_raw[i])) g_star[i] = 0.0;  // impute to mean
    else g_star[i] = (g_raw[i] - 2.0 * p) / sd;
  }
  return List::create(Named("g_star") = g_star, Named("maf") = p);
}

// OLS beta estimate and residuals
// [[Rcpp::export]]
List ols_resid(NumericVector r, NumericVector g_star) {
  int n = r.size();
  double ss_gg = 0.0, ss_rg = 0.0;
  for (int i = 0; i < n; i++) {
    ss_gg += g_star[i] * g_star[i];
    ss_rg += g_star[i] * r[i];
  }
  double beta = ss_rg / ss_gg;
  NumericVector e(n);
  for (int i = 0; i < n; i++) e[i] = r[i] - beta * g_star[i];
  return List::create(Named("e") = e, Named("beta") = beta);
}

// Apply RINT using precomputed phi_inv grid
// Returns RINT-transformed vector
// [[Rcpp::export]]
NumericVector apply_rint(NumericVector e, NumericVector phi_inv) {
  int n = e.size();
  // rank using stable sort
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
    [&e](int a, int b) { return e[a] < e[b]; });
  NumericVector e_rint(n);
  for (int i = 0; i < n; i++) e_rint[idx[i]] = phi_inv[i];
  return e_rint;
}

// Get group indices for K=3 genotype groups (round dosage to 0,1,2)
// [[Rcpp::export]]
List get_groups(NumericVector g_raw) {
  int n = g_raw.size();
  std::vector<std::vector<int>> grps(3);
  std::vector<int> n_g(3, 0);
  for (int i = 0; i < n; i++) {
    int g = (int)std::round(g_raw[i]);
    if (g < 0) g = 0;
    if (g > 2) g = 2;
    grps[g].push_back(i);
    n_g[g]++;
  }
  return List::create(
    Named("g0") = wrap(grps[0]),
    Named("g1") = wrap(grps[1]),
    Named("g2") = wrap(grps[2]),
    Named("n_g") = wrap(n_g)
  );
}

// ─────────────────────────────────────────────────────────────────────────────
// METHOD 1: LEVENE + RINT
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List levene_rint_cpp(NumericVector e_rint, IntegerVector grp_idx0,
                     IntegerVector grp_idx1, IntegerVector grp_idx2) {
  int n = e_rint.size();
  int K = 3;
  std::vector<IntegerVector> grps = {grp_idx0, grp_idx1, grp_idx2};
  std::vector<int> n_g = {(int)grp_idx0.size(),
                           (int)grp_idx1.size(),
                           (int)grp_idx2.size()};

  // Skip degenerate groups
  for (int k = 0; k < K; k++) if (n_g[k] < 2) {
    return List::create(Named("T_stat") = NA_REAL,
                        Named("pval")   = NA_REAL);
  }

  // Step 1: group means of e_rint
  std::vector<double> mu_g(K, 0.0);
  for (int k = 0; k < K; k++) {
    for (int i : grps[k]) mu_g[k] += e_rint[i];
    mu_g[k] /= n_g[k];
  }

  // Step 2: absolute deviations a_i
  std::vector<double> a(n);
  for (int k = 0; k < K; k++)
    for (int i : grps[k]) a[i] = std::abs(e_rint[i] - mu_g[k]);

  // Step 3: group means and grand mean of a
  std::vector<double> a_bar_g(K, 0.0);
  double a_bar = 0.0;
  for (int k = 0; k < K; k++) {
    for (int i : grps[k]) a_bar_g[k] += a[i];
    a_bar_g[k] /= n_g[k];
    a_bar += a_bar_g[k] * n_g[k];
  }
  a_bar /= n;

  // Step 4: F-statistic
  double SS_B = 0.0, SS_W = 0.0;
  for (int k = 0; k < K; k++) {
    SS_B += n_g[k] * std::pow(a_bar_g[k] - a_bar, 2.0);
    for (int i : grps[k]) SS_W += std::pow(a[i] - a_bar_g[k], 2.0);
  }
  if (SS_W < 1e-15) return List::create(Named("T_stat") = NA_REAL,
                                         Named("pval")   = NA_REAL);
  double F_stat = ((double)(n - K) / (double)(K - 1)) * SS_B / SS_W;
  double T_stat = (double)(K - 1) * F_stat;  // chi-sq approx with (K-1) df
  double pval   = R::pchisq(T_stat, K - 1, 0, 0);
  return List::create(Named("T_stat") = T_stat,
                      Named("pval")   = pval,
                      Named("F_stat") = F_stat);
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared helper: group median (used by BF and BF-noRINT)
// ─────────────────────────────────────────────────────────────────────────────

static double group_median(NumericVector v_in, IntegerVector idx) {
  int m = idx.size();
  std::vector<double> v(m);
  for (int j = 0; j < m; j++) v[j] = v_in[idx[j]];
  std::nth_element(v.begin(), v.begin() + m / 2, v.end());
  double med = v[m / 2];
  if (m % 2 == 0) {
    std::nth_element(v.begin(), v.begin() + m / 2 - 1, v.end());
    med = (med + v[m / 2 - 1]) / 2.0;
  }
  return med;
}

// ─────────────────────────────────────────────────────────────────────────────
// METHOD 1b: LEVENE WITHOUT RINT  (Wang et al. 2019 recommendation)
// Applied directly to OLS residuals; no rank transformation.
// Preserves Levene's original form but removes the RINT-induced
// contamination coupling.  Still susceptible to non-normality and
// the residual (non-RINT) contamination NCP at large n.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List levene_norint_cpp(NumericVector e_ols, IntegerVector grp_idx0,
                       IntegerVector grp_idx1, IntegerVector grp_idx2) {
  int n = e_ols.size();
  int K = 3;
  std::vector<IntegerVector> grps = {grp_idx0, grp_idx1, grp_idx2};
  std::vector<int> n_g = {(int)grp_idx0.size(),
                           (int)grp_idx1.size(),
                           (int)grp_idx2.size()};

  for (int k = 0; k < K; k++) if (n_g[k] < 2) {
    return List::create(Named("T_stat") = NA_REAL, Named("pval") = NA_REAL);
  }

  // Group means of OLS residuals
  std::vector<double> mu_g(K, 0.0);
  for (int k = 0; k < K; k++) {
    for (int i : grps[k]) mu_g[k] += e_ols[i];
    mu_g[k] /= n_g[k];
  }

  // Absolute deviations from group mean
  std::vector<double> a(n);
  for (int k = 0; k < K; k++)
    for (int i : grps[k]) a[i] = std::abs(e_ols[i] - mu_g[k]);

  // Group means and grand mean of absolute deviations
  std::vector<double> a_bar_g(K, 0.0);
  double a_bar = 0.0;
  for (int k = 0; k < K; k++) {
    for (int i : grps[k]) a_bar_g[k] += a[i];
    a_bar_g[k] /= n_g[k];
    a_bar += a_bar_g[k] * n_g[k];
  }
  a_bar /= n;

  double SS_B = 0.0, SS_W = 0.0;
  for (int k = 0; k < K; k++) {
    SS_B += n_g[k] * std::pow(a_bar_g[k] - a_bar, 2.0);
    for (int i : grps[k]) SS_W += std::pow(a[i] - a_bar_g[k], 2.0);
  }
  if (SS_W < 1e-15) return List::create(Named("T_stat") = NA_REAL,
                                         Named("pval")   = NA_REAL);
  double F_stat = ((double)(n - K) / (double)(K - 1)) * SS_B / SS_W;
  double T_stat = (double)(K - 1) * F_stat;
  double pval   = R::pchisq(T_stat, K - 1, 0, 0);
  return List::create(Named("T_stat") = T_stat,
                      Named("pval")   = pval,
                      Named("F_stat") = F_stat);
}

// ─────────────────────────────────────────────────────────────────────────────
// METHOD 2b: BROWN-FORSYTHE WITHOUT RINT  (median-based, no rank transform)
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List bf_norint_cpp(NumericVector e_ols, IntegerVector grp_idx0,
                   IntegerVector grp_idx1, IntegerVector grp_idx2) {
  int n = e_ols.size();
  int K = 3;
  std::vector<IntegerVector> grps = {grp_idx0, grp_idx1, grp_idx2};
  std::vector<int> n_g = {(int)grp_idx0.size(),
                           (int)grp_idx1.size(),
                           (int)grp_idx2.size()};

  for (int k = 0; k < K; k++) if (n_g[k] < 2) {
    return List::create(Named("T_stat") = NA_REAL, Named("pval") = NA_REAL);
  }

  // Group medians of OLS residuals
  std::vector<double> med_g(K);
  for (int k = 0; k < K; k++) med_g[k] = group_median(e_ols, grps[k]);

  // Absolute deviations from group median
  std::vector<double> a(n);
  for (int k = 0; k < K; k++)
    for (int i : grps[k]) a[i] = std::abs(e_ols[i] - med_g[k]);

  std::vector<double> a_bar_g(K, 0.0);
  double a_bar = 0.0;
  for (int k = 0; k < K; k++) {
    for (int i : grps[k]) a_bar_g[k] += a[i];
    a_bar_g[k] /= n_g[k];
    a_bar += a_bar_g[k] * n_g[k];
  }
  a_bar /= n;

  double SS_B = 0.0, SS_W = 0.0;
  for (int k = 0; k < K; k++) {
    SS_B += n_g[k] * std::pow(a_bar_g[k] - a_bar, 2.0);
    for (int i : grps[k]) SS_W += std::pow(a[i] - a_bar_g[k], 2.0);
  }
  if (SS_W < 1e-15) return List::create(Named("T_stat") = NA_REAL,
                                         Named("pval")   = NA_REAL);
  double F_stat = ((double)(n - K) / (double)(K - 1)) * SS_B / SS_W;
  double T_stat = (double)(K - 1) * F_stat;
  double pval   = R::pchisq(T_stat, K - 1, 0, 0);
  return List::create(Named("T_stat") = T_stat,
                      Named("pval")   = pval,
                      Named("F_stat") = F_stat);
}

// ─────────────────────────────────────────────────────────────────────────────
// ─────────────────────────────────────────────────────────────────────────────
// METHOD 2: BROWN-FORSYTHE + RINT
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List bf_rint_cpp(NumericVector e_rint, IntegerVector grp_idx0,
                 IntegerVector grp_idx1, IntegerVector grp_idx2) {
  int n = e_rint.size();
  int K = 3;
  std::vector<IntegerVector> grps = {grp_idx0, grp_idx1, grp_idx2};
  std::vector<int> n_g = {(int)grp_idx0.size(),
                           (int)grp_idx1.size(),
                           (int)grp_idx2.size()};

  for (int k = 0; k < K; k++) if (n_g[k] < 2) {
    return List::create(Named("T_stat") = NA_REAL, Named("pval") = NA_REAL);
  }

  // Group medians (BF differs from Levene only here)
  std::vector<double> med_g(K);
  for (int k = 0; k < K; k++) med_g[k] = group_median(e_rint, grps[k]);

  // Absolute deviations from median
  std::vector<double> a(n);
  for (int k = 0; k < K; k++)
    for (int i : grps[k]) a[i] = std::abs(e_rint[i] - med_g[k]);

  // Group means of a and grand mean
  std::vector<double> a_bar_g(K, 0.0);
  double a_bar = 0.0;
  for (int k = 0; k < K; k++) {
    for (int i : grps[k]) a_bar_g[k] += a[i];
    a_bar_g[k] /= n_g[k];
    a_bar += a_bar_g[k] * n_g[k];
  }
  a_bar /= n;

  double SS_B = 0.0, SS_W = 0.0;
  for (int k = 0; k < K; k++) {
    SS_B += n_g[k] * std::pow(a_bar_g[k] - a_bar, 2.0);
    for (int i : grps[k]) SS_W += std::pow(a[i] - a_bar_g[k], 2.0);
  }
  if (SS_W < 1e-15) return List::create(Named("T_stat") = NA_REAL,
                                         Named("pval")   = NA_REAL);
  double F_stat = ((double)(n - K) / (double)(K - 1)) * SS_B / SS_W;
  double T_stat = (double)(K - 1) * F_stat;
  double pval   = R::pchisq(T_stat, K - 1, 0, 0);
  return List::create(Named("T_stat") = T_stat, Named("pval") = pval);
}

// ─────────────────────────────────────────────────────────────────────────────
// METHOD 3: BARTLETT + RINT
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List bartlett_rint_cpp(NumericVector e_rint, IntegerVector grp_idx0,
                       IntegerVector grp_idx1, IntegerVector grp_idx2) {
  int n = e_rint.size();
  int K = 3;
  std::vector<IntegerVector> grps = {grp_idx0, grp_idx1, grp_idx2};
  std::vector<int> n_g = {(int)grp_idx0.size(),
                           (int)grp_idx1.size(),
                           (int)grp_idx2.size()};

  for (int k = 0; k < K; k++) if (n_g[k] < 2) {
    return List::create(Named("T_stat") = NA_REAL, Named("pval") = NA_REAL);
  }

  // One pass: group means and sums of squares
  std::vector<double> mu_g(K, 0.0), ss_g(K, 0.0);
  for (int k = 0; k < K; k++) {
    for (int i : grps[k]) mu_g[k] += e_rint[i];
    mu_g[k] /= n_g[k];
    for (int i : grps[k]) ss_g[k] += std::pow(e_rint[i] - mu_g[k], 2.0);
  }

  double dof_total = n - K;
  double sp2 = 0.0, sum_log_sg2 = 0.0;
  double c_inv = 0.0;
  for (int k = 0; k < K; k++) {
    double sg2 = ss_g[k] / (n_g[k] - 1.0);
    if (sg2 < 1e-15) return List::create(Named("T_stat") = NA_REAL,
                                          Named("pval")   = NA_REAL);
    sp2 += ss_g[k];
    sum_log_sg2 += (n_g[k] - 1.0) * std::log(sg2);
    c_inv += 1.0 / (n_g[k] - 1.0);
  }
  sp2 /= dof_total;
  double C_B = 1.0 + (c_inv - 1.0 / dof_total) / (3.0 * (K - 1.0));
  double B = (dof_total * std::log(sp2) - sum_log_sg2) / C_B;
  double pval = R::pchisq(B, K - 1, 0, 0);
  return List::create(Named("T_stat") = B, Named("pval") = pval);
}

// ─────────────────────────────────────────────────────────────────────────────
// METHOD 4: DGLM
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List dglm_cpp(NumericVector r, NumericVector g_star, int max_iter = 10) {
  int n = r.size();
  double ss_gg = 0.0;
  for (int i = 0; i < n; i++) ss_gg += g_star[i] * g_star[i];

  // Initial OLS
  double beta = 0.0;
  for (int i = 0; i < n; i++) beta += g_star[i] * r[i];
  beta /= ss_gg;

  std::vector<double> e(n), log_d2(n), log_s2(n), w(n);

  for (int iter = 0; iter < max_iter; iter++) {
    // Residuals and dispersion model
    double alpha = 0.0, ss_gg_w = 0.0;
    for (int i = 0; i < n; i++) {
      e[i] = r[i] - beta * g_star[i];
      log_d2[i] = std::log(std::max(e[i] * e[i], 1e-10));
    }
    for (int i = 0; i < n; i++) alpha += g_star[i] * log_d2[i];
    alpha /= ss_gg;
    for (int i = 0; i < n; i++) log_s2[i] = alpha * g_star[i];

    // Update mean model with precision weights
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; i++) {
      w[i] = std::exp(-log_s2[i]);
      num += w[i] * g_star[i] * r[i];
      den += w[i] * g_star[i] * g_star[i];
    }
    beta = (den > 1e-15) ? num / den : 0.0;
  }

  // Final residuals and score test
  double log_mean = 0.0;
  for (int i = 0; i < n; i++) {
    e[i] = r[i] - beta * g_star[i];
    log_d2[i] = std::log(std::max(e[i] * e[i], 1e-10));
    log_mean += log_d2[i];
  }
  log_mean /= n;

  double U = 0.0, I_info = 0.0;
  for (int i = 0; i < n; i++) {
    U += g_star[i] * (log_d2[i] - log_mean);
    I_info += g_star[i] * g_star[i];
  }
  I_info *= (M_PI * M_PI / 2.0);
  double T_D = (I_info > 1e-15) ? U * U / I_info : 0.0;
  double pval = R::pchisq(T_D, 1, 0, 0);
  return List::create(Named("T_stat") = T_D,
                      Named("pval")   = pval,
                      Named("beta")   = beta);
}

// ─────────────────────────────────────────────────────────────────────────────
// METHOD 5: POLYVAR (primary)
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List polyvar_scores_cpp(NumericVector e_rint, NumericVector g_star,
                        double beta_hat) {
  int n = e_rint.size();

  // Pass 1: absolute deviations and score accumulators
  double a_bar = 0.0;
  std::vector<double> a(n);
  double U1 = 0.0, U2 = 0.0, Umu = 0.0;
  for (int i = 0; i < n; i++) {
    a[i]   = std::abs(e_rint[i]);
    a_bar += a[i];
    Umu   += g_star[i] * e_rint[i];
  }
  a_bar /= n;

  // Pass 2: centred scores
  std::vector<double> u1v(n), u2v(n), umv(n);
  for (int i = 0; i < n; i++) {
    double ac  = a[i] - a_bar;
    double gs2 = g_star[i] * g_star[i];
    u1v[i] = g_star[i] * ac;
    u2v[i] = (gs2 - 1.0) * ac;
    umv[i] = g_star[i] * e_rint[i];
    U1    += u1v[i];
    U2    += u2v[i];
  }

  // Pass 3: empirical covariances for projection
  double cov1u = 0.0, cov2u = 0.0, var_u = 0.0;
  double sum_u1 = U1, sum_u2 = U2, sum_um = Umu;
  for (int i = 0; i < n; i++) {
    cov1u += u1v[i] * umv[i];
    cov2u += u2v[i] * umv[i];
    var_u += umv[i] * umv[i];
  }
  cov1u -= sum_u1 * sum_um / n;
  cov2u -= sum_u2 * sum_um / n;
  var_u -= sum_um * sum_um / n;

  if (std::abs(var_u) < 1e-15) {
    // No mean signal -- projection is identity
    // Fall back to unprojected scores
  }
  double c1 = (std::abs(var_u) > 1e-15) ? cov1u / var_u : 0.0;
  double c2 = (std::abs(var_u) > 1e-15) ? cov2u / var_u : 0.0;

  double U1c = U1 - c1 * Umu;
  double U2c = U2 - c2 * Umu;

  // Pass 4: variance-covariance of projected scores
  double V11 = 0.0, V22 = 0.0, V12 = 0.0;
  for (int i = 0; i < n; i++) {
    double u1ci = u1v[i] - c1 * umv[i];
    double u2ci = u2v[i] - c2 * umv[i];
    V11 += u1ci * u1ci;
    V22 += u2ci * u2ci;
    V12 += u1ci * u2ci;
  }

  double det = V11 * V22 - V12 * V12;
  if (det < 1e-15) return List::create(Named("T_C")      = NA_REAL,
                                        Named("pval")     = NA_REAL,
                                        Named("beta_hat") = beta_hat);

  double T_C = (V22 * U1c * U1c - 2.0 * V12 * U1c * U2c + V11 * U2c * U2c)
               / det;
  double pval = R::pchisq(T_C, 2, 0, 0);

  return List::create(Named("T_C")      = T_C,
                      Named("pval")     = pval,
                      Named("U1c")      = U1c,
                      Named("U2c")      = U2c,
                      Named("beta_hat") = beta_hat);
}

// LUT calibration: given T_C and beta_hat, return calibrated statistic
// lut_beta: sorted beta_abs grid; lut_q95: corresponding null q95 values
// [[Rcpp::export]]
List polyvar_calibrate_cpp(double T_C, double beta_hat,
                           NumericVector lut_beta,
                           NumericVector lut_q95) {
  // Linear interpolation of LUT
  int m = lut_beta.size();
  double ab = std::abs(beta_hat);
  ab = std::max(lut_beta[0], std::min(lut_beta[m-1], ab));

  // Find bracket
  int lo = 0;
  for (int k = 0; k < m - 1; k++) if (lut_beta[k] <= ab) lo = k;
  double t = (lut_beta[lo+1] > lut_beta[lo]) ?
    (ab - lut_beta[lo]) / (lut_beta[lo+1] - lut_beta[lo]) : 0.0;
  double q95 = lut_q95[lo] + t * (lut_q95[lo+1] - lut_q95[lo]);

  double chi2_q95 = R::qchisq(0.95, 2, 1, 0);  // = 5.991465
  double T_cal = (q95 > 1e-10) ? T_C * chi2_q95 / q95 : T_C;
  double pval  = R::pchisq(T_cal, 2, 0, 0);

  return List::create(Named("T_cal") = T_cal, Named("pval") = pval,
                      Named("q95_used") = q95);
}

// ─────────────────────────────────────────────────────────────────────────────
// METHOD 6: POLYVAR-LITE
// ─────────────────────────────────────────────────────────────────────────────

// Compute bin boundaries from LOCO residuals (one-time per phenotype)
// [[Rcpp::export]]
List polyvar_lite_setup_cpp(NumericVector r, int Q = 100) {
  int n = r.size();
  // Compute quantile breaks
  NumericVector probs(Q + 1);
  for (int j = 0; j <= Q; j++) probs[j] = (double)j / Q;
  // Sort r to get empirical quantiles
  std::vector<double> rs(r.begin(), r.end());
  std::sort(rs.begin(), rs.end());

  NumericVector breaks(Q + 1);
  breaks[0] = R_NegInf;
  breaks[Q] = R_PosInf;
  for (int j = 1; j < Q; j++) {
    double p = probs[j];
    double pos = p * (n - 1);
    int lo = (int)pos;
    double frac = pos - lo;
    breaks[j] = rs[lo] + frac * (rs[std::min(lo+1, n-1)] - rs[lo]);
  }

  // Bin-centre normal quantiles
  NumericVector bin_q(Q);
  for (int j = 0; j < Q; j++) {
    double pmid = (probs[j] + probs[j + 1]) / 2.0;
    pmid = std::max(1e-8, std::min(1.0 - 1e-8, pmid));
    bin_q[j] = R::qnorm(pmid, 0.0, 1.0, 1, 0);
  }

  return List::create(Named("breaks") = breaks, Named("bin_q") = bin_q,
                      Named("Q") = Q);
}

// Per-SNP PolyVAR-lite: bin assignment + same scoring as polyvar
// [[Rcpp::export]]
List polyvar_lite_stat_cpp(NumericVector r, NumericVector g_star,
                           double beta_hat,
                           NumericVector breaks, NumericVector bin_q) {
  int n = r.size();
  int Q = bin_q.size();

  // OLS residuals
  std::vector<double> e(n);
  for (int i = 0; i < n; i++) e[i] = r[i] - beta_hat * g_star[i];

  // Bin assignment via binary search: O(n log Q)
  NumericVector e_rint(n);
  for (int i = 0; i < n; i++) {
    int lo = 0, hi = Q - 1;
    while (lo < hi) {
      int mid = (lo + hi + 1) / 2;
      if (e[i] >= breaks[mid]) lo = mid;
      else hi = mid - 1;
    }
    e_rint[i] = bin_q[lo];
  }

  // Same scoring as polyvar_scores_cpp
  return polyvar_scores_cpp(e_rint, g_star, beta_hat);
}

// ─────────────────────────────────────────────────────────────────────────────
// BATCH PROCESSING: all methods, multiple SNPs
// Takes matrix G (n x M) and vector r (n), phi_inv, LUT
// Returns data.frame with T_stat, pval, beta, gamma2 per SNP
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
DataFrame vqtl_batch_cpp(NumericMatrix G,        // n x M genotype matrix
                         NumericVector r,         // n phenotype residuals
                         NumericVector phi_inv,   // n RINT grid
                         std::string method,      // method name
                         NumericVector lut_beta,  // LUT for PolyVAR
                         NumericVector lut_q95,   // LUT for PolyVAR
                         int Q_lite = 100,        // PolyVAR-lite bins
                         int dglm_iter = 10) {    // DGLM iterations

  int n = G.nrow();
  int M = G.ncol();

  NumericVector T_stat(M, NA_REAL);
  NumericVector pval(M, NA_REAL);
  NumericVector beta_hat(M, NA_REAL);
  NumericVector gamma2_hat(M, NA_REAL);

  // Per-phenotype lite setup
  List lite_setup;
  if (method == "polyvar_lite") {
    lite_setup = polyvar_lite_setup_cpp(r, Q_lite);
  }

  for (int j = 0; j < M; j++) {
    NumericVector g_raw = G(_, j);

    // Standardise genotype
    List sg = standardise_geno(g_raw);
    NumericVector g_star = sg["g_star"];
    double maf = sg["maf"];

    // OLS mean removal
    List ol = ols_resid(r, g_star);
    NumericVector e = ol["e"];
    double bhat = ol["beta"];
    beta_hat[j] = bhat;

    // RINT
    NumericVector e_rint = apply_rint(e, phi_inv);

    // Group indices (for Levene/BF/Bartlett variants)
    List grp;
    if (method == "levene" || method == "bf" || method == "bartlett" ||
        method == "levene_norint" || method == "bf_norint") {
      grp = get_groups(g_raw);
    }

    List result;

    if (method == "levene") {
      result = levene_rint_cpp(e_rint,
        grp["g0"], grp["g1"], grp["g2"]);
    } else if (method == "levene_norint") {
      result = levene_norint_cpp(e,
        grp["g0"], grp["g1"], grp["g2"]);
    } else if (method == "bf") {
      result = bf_rint_cpp(e_rint,
        grp["g0"], grp["g1"], grp["g2"]);
    } else if (method == "bf_norint") {
      result = bf_norint_cpp(e,
        grp["g0"], grp["g1"], grp["g2"]);
    } else if (method == "bartlett") {
      result = bartlett_rint_cpp(e_rint,
        grp["g0"], grp["g1"], grp["g2"]);
    } else if (method == "dglm") {
      result = dglm_cpp(r, g_star, dglm_iter);
    } else if (method == "polyvar") {
      List sc = polyvar_scores_cpp(e_rint, g_star, bhat);
      if (!NumericVector::is_na(as<double>(sc["T_C"]))) {
        result = polyvar_calibrate_cpp(sc["T_C"], bhat,
                                       lut_beta, lut_q95);
        result["beta_hat"] = bhat;
      } else {
        result = List::create(Named("T_cal") = NA_REAL,
                              Named("pval")  = NA_REAL);
      }
    } else if (method == "polyvar_lite") {
      NumericVector breaks = lite_setup["breaks"];
      NumericVector bin_q  = lite_setup["bin_q"];
      List sc = polyvar_lite_stat_cpp(r, g_star, bhat, breaks, bin_q);
      if (!NumericVector::is_na(as<double>(sc["T_C"]))) {
        result = polyvar_calibrate_cpp(sc["T_C"], bhat,
                                       lut_beta, lut_q95);
        result["beta_hat"] = bhat;
      } else {
        result = List::create(Named("T_cal") = NA_REAL,
                              Named("pval")  = NA_REAL);
      }
    }

    // Extract stats
    double T_j = NA_REAL, p_j = NA_REAL;
    if (result.containsElementNamed("T_cal")) {
      T_j = as<double>(result["T_cal"]);
      p_j = as<double>(result["pval"]);
    } else if (result.containsElementNamed("T_stat")) {
      T_j = as<double>(result["T_stat"]);
      p_j = as<double>(result["pval"]);
    }

    T_stat[j] = T_j;
    pval[j]   = p_j;
    int df_j  = (method == "dglm") ? 1 : 2;
    if (!NumericVector::is_na(T_j) && T_j > 0) {
      double kap = 1.0 / (2.0 * maf * (1.0 - maf)) - 1.0;
      double phi = 0.438 * n * kap;
      gamma2_hat[j] = std::sqrt(std::max(T_j - df_j, 0.0) /
                                 std::max(phi, 1.0));
    }
  }

  return DataFrame::create(
    Named("T_stat")     = T_stat,
    Named("pval")       = pval,
    Named("beta_hat")   = beta_hat,
    Named("gamma2_hat") = gamma2_hat
  );
}

// polyvar_fast.cpp
// Optimised per-SNP PolyVAR and PolyVAR-lite fully in C++.
// Key design: one self-contained function per method, no R callbacks.
// OLS + sort/binning + score accumulation all in one pass where possible.
//
// Benchmark target: polyvar_full_cpp and polyvar_lite_cpp should each
// beat their R-loop counterparts by the full sort-cost ratio.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>

using namespace Rcpp;

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers (not exported)
// ─────────────────────────────────────────────────────────────────────────────

// Score projection + Mahalanobis T_C from (U1,U2,Umu) and raw score vectors
static double compute_TC(double U1, double U2, double Umu,
                          const std::vector<double>& u1v,
                          const std::vector<double>& u2v,
                          const std::vector<double>& umv, int n) {
  double cov1u = 0, cov2u = 0, var_u = 0;
  double s1 = U1, s2 = U2, sm = Umu;
  for (int i = 0; i < n; i++) {
    cov1u += u1v[i] * umv[i];
    cov2u += u2v[i] * umv[i];
    var_u += umv[i] * umv[i];
  }
  cov1u -= s1 * sm / n;
  cov2u -= s2 * sm / n;
  var_u -= sm * sm / n;

  double c1 = (var_u > 1e-15) ? cov1u / var_u : 0.0;
  double c2 = (var_u > 1e-15) ? cov2u / var_u : 0.0;
  double U1c = U1 - c1 * Umu;
  double U2c = U2 - c2 * Umu;

  double V11 = 0, V22 = 0, V12 = 0;
  for (int i = 0; i < n; i++) {
    double v1 = u1v[i] - c1 * umv[i];
    double v2 = u2v[i] - c2 * umv[i];
    V11 += v1 * v1;  V22 += v2 * v2;  V12 += v1 * v2;
  }
  double det = V11 * V22 - V12 * V12;
  if (det < 1e-15) return R_NaReal;
  return (V22*U1c*U1c - 2.0*V12*U1c*U2c + V11*U2c*U2c) / det;
}

// ─────────────────────────────────────────────────────────────────────────────
// POLYVAR FULL (exact RINT via C++ std::sort)
// Self-contained: genotype + residuals -> T_C
// No R callbacks after entry.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List polyvar_full_cpp(NumericVector r,       // n phenotype residuals
                      NumericVector g_raw,   // n raw dosages {0,1,2}
                      NumericVector phi_inv) // n RINT grid (precomputed once)
{
  int n = r.size();

  // Step 1: standardise genotype                               O(n)
  double sum_g = 0; int n_valid = 0;
  for (int i = 0; i < n; i++)
    if (!NumericVector::is_na(g_raw[i])) { sum_g += g_raw[i]; n_valid++; }
  double p   = std::max(1e-6, std::min(1-1e-6, sum_g / (2.0*n_valid)));
  double psd = std::sqrt(2.0*p*(1.0-p));
  std::vector<double> gstar(n);
  for (int i = 0; i < n; i++)
    gstar[i] = NumericVector::is_na(g_raw[i]) ? 0.0 : (g_raw[i]-2*p)/psd;

  // Step 2: OLS beta + residual                                O(n)
  double sg2=0, srg=0;
  for (int i=0;i<n;i++){sg2+=gstar[i]*gstar[i]; srg+=gstar[i]*r[i];}
  double beta = srg/sg2;
  std::vector<double> e(n);
  for (int i=0;i<n;i++) e[i]=r[i]-beta*gstar[i];

  // Step 3: RINT via C++ sort                                  O(n log n)
  std::vector<int> idx(n);
  std::iota(idx.begin(),idx.end(),0);
  std::stable_sort(idx.begin(),idx.end(),[&e](int a,int b){return e[a]<e[b];});
  std::vector<double> er(n);
  for (int i=0;i<n;i++) er[idx[i]] = phi_inv[i];

  // Step 4: score accumulation                                 O(n)
  double a_bar=0, U1=0, U2=0, Umu=0;
  std::vector<double> a(n);
  for (int i=0;i<n;i++){a[i]=std::abs(er[i]); a_bar+=a[i]; Umu+=gstar[i]*er[i];}
  a_bar /= n;
  std::vector<double> u1v(n),u2v(n),umv(n);
  for (int i=0;i<n;i++){
    double ac=a[i]-a_bar; double gs2=gstar[i]*gstar[i];
    u1v[i]=gstar[i]*ac;  u2v[i]=(gs2-1)*ac;  umv[i]=gstar[i]*er[i];
    U1+=u1v[i]; U2+=u2v[i];
  }

  double TC = compute_TC(U1,U2,Umu,u1v,u2v,umv,n);
  double pv = NumericVector::is_na(TC) ? NA_REAL : R::pchisq(TC,2,0,0);
  return List::create(Named("T_C")=TC, Named("pval")=pv,
                      Named("beta")=beta, Named("maf")=p);
}

// ─────────────────────────────────────────────────────────────────────────────
// POLYVAR-LITE SETUP: bin boundaries from LOCO residuals (C++, O(n log n))
// Called ONCE per phenotype; not per SNP.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List polyvar_lite_setup_fast(NumericVector r, int Q = 100) {
  int n = r.size();
  // Sort once to get empirical quantile breakpoints              O(n log n)
  std::vector<double> rs(r.begin(), r.end());
  std::sort(rs.begin(), rs.end());

  NumericVector breaks(Q+1), bin_q(Q);
  breaks[0] = R_NegInf; breaks[Q] = R_PosInf;
  for (int j=1; j<Q; j++) {
    double pos = (double)j/Q * (n-1);
    int lo = (int)pos; double frac = pos-lo;
    breaks[j] = rs[lo] + frac*(rs[std::min(lo+1,n-1)]-rs[lo]);
  }
  for (int j=0; j<Q; j++) {
    double pmid = ((double)j+0.5)/Q;
    pmid = std::max(1e-8, std::min(1-1e-8, pmid));
    bin_q[j] = R::qnorm(pmid,0,1,1,0);
  }
  return List::create(Named("breaks")=breaks, Named("bin_q")=bin_q, Named("Q")=Q);
}

// ─────────────────────────────────────────────────────────────────────────────
// POLYVAR-LITE FAST: per-SNP O(n log Q) — no sort, only binary search
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List polyvar_lite_fast(NumericVector r,        // n phenotype residuals
                       NumericVector g_raw,    // n raw dosages
                       NumericVector breaks,   // Q+1 bin boundaries
                       NumericVector bin_q)    // Q bin-centre quantiles
{
  int n = r.size();
  int Q = bin_q.size();

  // Step 1: standardise genotype                               O(n)
  double sum_g=0; int n_valid=0;
  for (int i=0;i<n;i++)
    if (!NumericVector::is_na(g_raw[i])){sum_g+=g_raw[i]; n_valid++;}
  double p   = std::max(1e-6,std::min(1-1e-6,sum_g/(2.0*n_valid)));
  double psd = std::sqrt(2.0*p*(1.0-p));
  std::vector<double> gstar(n);
  for (int i=0;i<n;i++)
    gstar[i] = NumericVector::is_na(g_raw[i]) ? 0.0 : (g_raw[i]-2*p)/psd;

  // Step 2: OLS                                                O(n)
  double sg2=0,srg=0;
  for (int i=0;i<n;i++){sg2+=gstar[i]*gstar[i]; srg+=gstar[i]*r[i];}
  double beta=srg/sg2;

  // Step 3: bin assignment via binary search on e              O(n log Q)
  // Simultaneously accumulate scores in the same pass
  double a_bar=0, U1=0, U2=0, Umu=0;
  std::vector<double> a(n), er(n), u1v(n), u2v(n), umv(n);

  for (int i=0;i<n;i++){
    double ei = r[i]-beta*gstar[i];
    // Binary search: find j such that breaks[j] <= ei < breaks[j+1]
    int lo=0, hi=Q-1;
    while (lo<hi){
      int mid=(lo+hi+1)/2;
      if (ei>=(double)breaks[mid]) lo=mid; else hi=mid-1;
    }
    er[i] = bin_q[lo];
    a[i]  = std::abs(er[i]);
    a_bar += a[i];
    Umu   += gstar[i]*er[i];
  }
  a_bar /= n;

  // Score accumulation pass                                    O(n)
  for (int i=0;i<n;i++){
    double ac=a[i]-a_bar; double gs2=gstar[i]*gstar[i];
    u1v[i]=gstar[i]*ac; u2v[i]=(gs2-1)*ac; umv[i]=gstar[i]*er[i];
    U1+=u1v[i]; U2+=u2v[i];
  }

  double TC = compute_TC(U1,U2,Umu,u1v,u2v,umv,n);
  double pv = NumericVector::is_na(TC) ? NA_REAL : R::pchisq(TC,2,0,0);
  return List::create(Named("T_C")=TC, Named("pval")=pv,
                      Named("beta")=beta, Named("maf")=p);
}

// ─────────────────────────────────────────────────────────────────────────────
// BATCH: genome-wide scan, all SNPs for one phenotype, fully in C++
// Returns data.frame. Accepts phi_inv for exact PolyVAR or
// breaks/bin_q for PolyVAR-lite.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
DataFrame polyvar_scan_cpp(NumericMatrix G,        // n x M dosage matrix
                            NumericVector r,         // n phenotype residuals
                            NumericVector phi_inv,   // n RINT grid (PolyVAR exact)
                            NumericVector breaks,    // Q+1 bins (PolyVAR-lite)
                            NumericVector bin_q,     // Q bin centres
                            bool use_lite,           // true = PolyVAR-lite
                            NumericVector lut_beta,  // LUT x-axis
                            NumericVector lut_q95,   // LUT y-axis
                            double min_maf = 0.01) {
  int n = G.nrow();
  int M = G.ncol();
  int Q = bin_q.size();
  int lut_m = lut_beta.size();
  double chi2_q95 = R::qchisq(0.95, 2, 1, 0);

  NumericVector Tv(M,NA_REAL), pv(M,NA_REAL), bv(M,NA_REAL), gv(M,NA_REAL);
  NumericVector maf_v(M,NA_REAL);

  for (int j=0; j<M; j++) {
    // Step 1: standardise genotype
    double sum_g=0; int n_valid=0;
    for (int i=0;i<n;i++){
      double gi=G(i,j);
      if (!ISNA(gi)){sum_g+=gi; n_valid++;}
    }
    double p = sum_g/(2.0*n_valid);
    if (p<min_maf || p>(1-min_maf)) continue;
    double psd=std::sqrt(2*p*(1-p));
    double kap=1.0/(2*p*(1-p))-1.0;

    std::vector<double> gstar(n);
    for (int i=0;i<n;i++) gstar[i]=ISNA(G(i,j))?0.0:(G(i,j)-2*p)/psd;

    // Step 2: OLS
    double sg2=0,srg=0;
    for (int i=0;i<n;i++){sg2+=gstar[i]*gstar[i]; srg+=gstar[i]*r[i];}
    double beta=srg/sg2;

    // Step 3 + 4: RINT or bin assignment + score accumulation
    std::vector<double> er(n), a(n), u1v(n), u2v(n), umv(n);
    double a_bar=0, U1=0, U2=0, Umu=0;

    if (!use_lite) {
      // Exact RINT via C++ sort                                O(n log n)
      std::vector<int> idx(n); std::iota(idx.begin(),idx.end(),0);
      std::vector<double> e(n);
      for (int i=0;i<n;i++) e[i]=r[i]-beta*gstar[i];
      std::stable_sort(idx.begin(),idx.end(),[&e](int a,int b){return e[a]<e[b];});
      for (int i=0;i<n;i++) er[idx[i]]=phi_inv[i];
    } else {
      // Bin assignment O(n log Q)
      for (int i=0;i<n;i++){
        double ei=r[i]-beta*gstar[i];
        int lo=0,hi=Q-1;
        while(lo<hi){int mid=(lo+hi+1)/2; if(ei>=(double)breaks[mid])lo=mid;else hi=mid-1;}
        er[i]=bin_q[lo];
      }
    }

    // Score accumulation                                       O(n)
    for (int i=0;i<n;i++){a[i]=std::abs(er[i]); a_bar+=a[i]; Umu+=gstar[i]*er[i];}
    a_bar/=n;
    for (int i=0;i<n;i++){
      double ac=a[i]-a_bar; double gs2=gstar[i]*gstar[i];
      u1v[i]=gstar[i]*ac; u2v[i]=(gs2-1)*ac; umv[i]=gstar[i]*er[i];
      U1+=u1v[i]; U2+=u2v[i];
    }

    double TC = compute_TC(U1,U2,Umu,u1v,u2v,umv,n);
    if (ISNA(TC)) continue;

    // LUT calibration                                          O(log lut_m)
    double ab=std::abs(beta);
    ab=std::max((double)lut_beta[0],std::min((double)lut_beta[lut_m-1],ab));
    int lo=0;
    for (int k=0;k<lut_m-1;k++) if(lut_beta[k]<=ab) lo=k;
    double t=(lut_beta[lo+1]>lut_beta[lo])?(ab-lut_beta[lo])/(lut_beta[lo+1]-lut_beta[lo]):0.0;
    double q95=lut_q95[lo]+t*(lut_q95[lo+1]-lut_q95[lo]);
    double TC_cal=(q95>1e-10)?TC*chi2_q95/q95:TC;
    double pvj=R::pchisq(TC_cal,2,0,0);

    Tv[j]=TC_cal; pv[j]=pvj; bv[j]=beta; maf_v[j]=p;
    gv[j]=std::sqrt(std::max(TC_cal-2.0,0.0)/std::max(0.438*n*kap,1.0));
  }

  return DataFrame::create(Named("T_stat")=Tv, Named("pval")=pv,
                            Named("beta")=bv,   Named("gamma2")=gv,
                            Named("maf")=maf_v);
}

// ─────────────────────────────────────────────────────────────────────────────
// FAST BATCH for all methods (Levene, BF, Bartlett also in single-pass C++)
// ─────────────────────────────────────────────────────────────────────────────

// Internal: Levene F stat on pre-RINT residuals + group arrays
static double levene_F(const std::vector<double>& er,
                        const std::vector<int>* grp,
                        const std::vector<int>& ng, int n, int K) {
  std::vector<double> mu(K,0), ag(K,0); double a_grand=0;
  std::vector<double> a(n);
  for (int k=0;k<K;k++){
    for (int i:grp[k]) mu[k]+=er[i];
    mu[k]/=ng[k];
  }
  for (int k=0;k<K;k++) for (int i:grp[k]){a[i]=std::abs(er[i]-mu[k]); ag[k]+=a[i];}
  for (int k=0;k<K;k++){ag[k]/=ng[k]; a_grand+=ag[k]*ng[k];}
  a_grand/=n;
  double ssb=0,ssw=0;
  for (int k=0;k<K;k++){
    ssb+=ng[k]*std::pow(ag[k]-a_grand,2);
    for (int i:grp[k]) ssw+=std::pow(a[i]-ag[k],2);
  }
  if (ssw<1e-15) return 0;
  return ((double)(n-K)/(K-1))*ssb/ssw;
}

// [[Rcpp::export]]
DataFrame all_methods_scan_cpp(NumericMatrix G,
                                NumericVector r,
                                NumericVector phi_inv,
                                NumericVector breaks,
                                NumericVector bin_q,
                                NumericVector lut_beta,
                                NumericVector lut_q95,
                                double min_maf = 0.01) {
  int n=G.nrow(), M=G.ncol(), Q=bin_q.size(), lut_m=lut_beta.size();
  double chi2_q95=R::qchisq(0.95,2,1,0);
  int nmeth=6; // polyvar, polyvar_lite, levene, bf, bartlett, dglm
  // Store results: T and p per method per SNP
  NumericMatrix Tm(M,nmeth), Pm(M,nmeth), Bm(M,nmeth);
  for (int jj=0;jj<M*nmeth;jj++){Tm[jj]=NA_REAL;Pm[jj]=NA_REAL;}

  for (int j=0;j<M;j++){
    double sum_g=0; int n_valid=0;
    for (int i=0;i<n;i++){double gi=G(i,j); if(!ISNA(gi)){sum_g+=gi;n_valid++;}}
    double p=sum_g/(2.0*n_valid);
    if(p<min_maf||p>(1-min_maf)) continue;
    double psd=std::sqrt(2*p*(1-p));

    std::vector<double> gstar(n);
    for(int i=0;i<n;i++) gstar[i]=ISNA(G(i,j))?0.0:(G(i,j)-2*p)/psd;

    double sg2=0,srg=0;
    for(int i=0;i<n;i++){sg2+=gstar[i]*gstar[i]; srg+=gstar[i]*r[i];}
    double beta=srg/sg2;
    for(int i=0;i<n;i++) Bm(j,0)=beta; // store once

    // Group indices (for Levene/BF/Bartlett)
    std::vector<int> grp[3]; std::vector<int> ng(3,0);
    for(int i=0;i<n;i++){int g=std::round(G(i,j)); if(g<0)g=0;if(g>2)g=2; grp[g].push_back(i); ng[g]++;}
    bool grp_ok=true; for(int k=0;k<3;k++) if(ng[k]<2){grp_ok=false;break;}

    // Exact RINT (for polyvar, levene, bf, bartlett)
    std::vector<double> e(n), er(n);
    for(int i=0;i<n;i++) e[i]=r[i]-beta*gstar[i];
    std::vector<int> idx(n); std::iota(idx.begin(),idx.end(),0);
    std::stable_sort(idx.begin(),idx.end(),[&e](int a,int b){return e[a]<e[b];});
    for(int i=0;i<n;i++) er[idx[i]]=phi_inv[i];

    // Lite RINT
    std::vector<double> er_lite(n);
    for(int i=0;i<n;i++){
      double ei=r[i]-beta*gstar[i]; int lo=0,hi=Q-1;
      while(lo<hi){int mid=(lo+hi+1)/2;if(ei>=(double)breaks[mid])lo=mid;else hi=mid-1;}
      er_lite[i]=bin_q[lo];
    }

    // ── PolyVAR (exact) ──────────────────────────────────────────
    {
      std::vector<double> a(n),u1v(n),u2v(n),umv(n);
      double ab=0,U1=0,U2=0,Um=0;
      for(int i=0;i<n;i++){a[i]=std::abs(er[i]);ab+=a[i];Um+=gstar[i]*er[i];}
      ab/=n;
      for(int i=0;i<n;i++){
        double ac=a[i]-ab,gs2=gstar[i]*gstar[i];
        u1v[i]=gstar[i]*ac; u2v[i]=(gs2-1)*ac; umv[i]=gstar[i]*er[i];
        U1+=u1v[i]; U2+=u2v[i];
      }
      double TC=compute_TC(U1,U2,Um,u1v,u2v,umv,n);
      if(!ISNA(TC)){
        double ab2=std::abs(beta); ab2=std::max((double)lut_beta[0],std::min((double)lut_beta[lut_m-1],ab2));
        int lo=0; for(int k=0;k<lut_m-1;k++) if(lut_beta[k]<=ab2) lo=k;
        double t2=(lut_beta[lo+1]>lut_beta[lo])?(ab2-lut_beta[lo])/(lut_beta[lo+1]-lut_beta[lo]):0;
        double q95=lut_q95[lo]+t2*(lut_q95[lo+1]-lut_q95[lo]);
        double TC2=(q95>1e-10)?TC*chi2_q95/q95:TC;
        Tm(j,0)=TC2; Pm(j,0)=R::pchisq(TC2,2,0,0);
      }
    }

    // ── PolyVAR-lite ─────────────────────────────────────────────
    {
      std::vector<double> a(n),u1v(n),u2v(n),umv(n);
      double ab=0,U1=0,U2=0,Um=0;
      for(int i=0;i<n;i++){a[i]=std::abs(er_lite[i]);ab+=a[i];Um+=gstar[i]*er_lite[i];}
      ab/=n;
      for(int i=0;i<n;i++){
        double ac=a[i]-ab,gs2=gstar[i]*gstar[i];
        u1v[i]=gstar[i]*ac; u2v[i]=(gs2-1)*ac; umv[i]=gstar[i]*er_lite[i];
        U1+=u1v[i]; U2+=u2v[i];
      }
      double TC=compute_TC(U1,U2,Um,u1v,u2v,umv,n);
      if(!ISNA(TC)){
        double ab2=std::abs(beta); ab2=std::max((double)lut_beta[0],std::min((double)lut_beta[lut_m-1],ab2));
        int lo=0; for(int k=0;k<lut_m-1;k++) if(lut_beta[k]<=ab2) lo=k;
        double t2=(lut_beta[lo+1]>lut_beta[lo])?(ab2-lut_beta[lo])/(lut_beta[lo+1]-lut_beta[lo]):0;
        double q95=lut_q95[lo]+t2*(lut_q95[lo+1]-lut_q95[lo]);
        double TC2=(q95>1e-10)?TC*chi2_q95/q95:TC;
        Tm(j,1)=TC2; Pm(j,1)=R::pchisq(TC2,2,0,0);
      }
    }

    if(!grp_ok) continue;

    // ── Levene ───────────────────────────────────────────────────
    {
      double F=levene_F(er,grp,ng,n,3);
      double T=(3-1)*F; Tm(j,2)=T; Pm(j,2)=R::pchisq(T,2,0,0);
    }

    // ── BF (median-based) ────────────────────────────────────────
    {
      std::vector<double> med(3); std::vector<double> a(n);
      for(int k=0;k<3;k++){
        int m=ng[k]; std::vector<double> v(m);
        for(int jj=0;jj<m;jj++) v[jj]=er[grp[k][jj]];
        std::nth_element(v.begin(),v.begin()+m/2,v.end());
        med[k]=v[m/2];
        if(m%2==0){std::nth_element(v.begin(),v.begin()+m/2-1,v.end()); med[k]=(med[k]+v[m/2-1])/2;}
      }
      for(int k=0;k<3;k++) for(int i:grp[k]) a[i]=std::abs(er[i]-med[k]);
      std::vector<double> ag(3,0); double agrand=0;
      for(int k=0;k<3;k++){for(int i:grp[k])ag[k]+=a[i]; ag[k]/=ng[k]; agrand+=ag[k]*ng[k];}
      agrand/=n;
      double ssb=0,ssw=0;
      for(int k=0;k<3;k++){ssb+=ng[k]*std::pow(ag[k]-agrand,2); for(int i:grp[k])ssw+=std::pow(a[i]-ag[k],2);}
      double F=(ssw>1e-15)?((double)(n-3)/2)*ssb/ssw:0;
      double T=2*F; Tm(j,3)=T; Pm(j,3)=R::pchisq(T,2,0,0);
    }

    // ── Bartlett ─────────────────────────────────────────────────
    {
      std::vector<double> mu(3,0),ss(3,0);
      for(int k=0;k<3;k++){for(int i:grp[k])mu[k]+=er[i]; mu[k]/=ng[k]; for(int i:grp[k])ss[k]+=std::pow(er[i]-mu[k],2);}
      double sp2=0,slogs=0,cinv=0;
      for(int k=0;k<3;k++){
        double sg2=ss[k]/(ng[k]-1); if(sg2<1e-15) goto bartlett_skip;
        sp2+=ss[k]; slogs+=(ng[k]-1)*std::log(sg2); cinv+=1.0/(ng[k]-1);
      }
      sp2/=(n-3);
      {double CB=1+(cinv-1.0/(n-3))/(3*2.0); double B=(((n-3)*std::log(sp2)-slogs)/CB); Tm(j,4)=B; Pm(j,4)=R::pchisq(B,2,0,0);}
      goto bartlett_done;
      bartlett_skip:; bartlett_done:;
    }

    // ── DGLM ─────────────────────────────────────────────────────
    {
      double b2=beta;
      std::vector<double> ei(n),ld(n),ls(n),w(n);
      for(int iter=0;iter<10;iter++){
        double al=0,ss2=0;
        for(int i=0;i<n;i++){ei[i]=r[i]-b2*gstar[i]; ld[i]=std::log(std::max(ei[i]*ei[i],1e-10));}
        for(int i=0;i<n;i++){al+=gstar[i]*ld[i]; ss2+=gstar[i]*gstar[i];}
        al/=ss2; for(int i=0;i<n;i++) ls[i]=al*gstar[i];
        double num=0,den=0;
        for(int i=0;i<n;i++){w[i]=std::exp(-ls[i]); num+=w[i]*gstar[i]*r[i]; den+=w[i]*gstar[i]*gstar[i];}
        b2=(den>1e-15)?num/den:0;
      }
      double lmean=0,U=0,II=0;
      for(int i=0;i<n;i++){ei[i]=r[i]-b2*gstar[i]; ld[i]=std::log(std::max(ei[i]*ei[i],1e-10)); lmean+=ld[i];}
      lmean/=n;
      for(int i=0;i<n;i++){U+=gstar[i]*(ld[i]-lmean); II+=gstar[i]*gstar[i];}
      II*=(M_PI*M_PI/2);
      double TD=(II>1e-15)?U*U/II:0; Tm(j,5)=TD; Pm(j,5)=R::pchisq(TD,1,0,0);
    }
  }

  // Pack into named vectors
  NumericVector pv(M),pl(M),plev(M),pbf(M),pbar(M),pdg(M);
  NumericVector tv(M),tl(M),tlev(M),tbf(M),tbar(M),tdg(M);
  for(int j=0;j<M;j++){
    tv[j]=Tm(j,0); pv[j]=Pm(j,0);
    tl[j]=Tm(j,1); pl[j]=Pm(j,1);
    tlev[j]=Tm(j,2); plev[j]=Pm(j,2);
    tbf[j]=Tm(j,3); pbf[j]=Pm(j,3);
    tbar[j]=Tm(j,4); pbar[j]=Pm(j,4);
    tdg[j]=Tm(j,5); pdg[j]=Pm(j,5);
  }

  return DataFrame::create(
    Named("T_polyvar")=tv,      Named("p_polyvar")=pv,
    Named("T_lite")=tl,         Named("p_lite")=pl,
    Named("T_levene")=tlev,     Named("p_levene")=plev,
    Named("T_bf")=tbf,          Named("p_bf")=pbf,
    Named("T_bartlett")=tbar,   Named("p_bartlett")=pbar,
    Named("T_dglm")=tdg,        Named("p_dglm")=pdg
  );
}

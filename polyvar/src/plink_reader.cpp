// plink_reader.cpp
// Streaming reader for PLINK BED/BIM/FAM files
// Loads genotype dosages for a batch of SNPs
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <cmath>

using namespace Rcpp;

// Read the PLINK BIM file, return data frame of SNP info
// [[Rcpp::export]]
DataFrame read_bim_cpp(std::string bim_file) {
  std::ifstream f(bim_file);
  if (!f.is_open()) stop("Cannot open BIM file: " + bim_file);

  std::vector<std::string> chr_v, snpid_v, a1_v, a2_v;
  std::vector<double> pos_v, cm_v;

  std::string chr, snpid, a1, a2;
  double cm, pos;
  while (f >> chr >> snpid >> cm >> pos >> a1 >> a2) {
    chr_v.push_back(chr);  snpid_v.push_back(snpid);
    cm_v.push_back(cm);    pos_v.push_back(pos);
    a1_v.push_back(a1);    a2_v.push_back(a2);
  }
  f.close();
  return DataFrame::create(
    Named("chr")  = chr_v,
    Named("snpid")= snpid_v,
    Named("cm")   = cm_v,
    Named("pos")  = pos_v,
    Named("a1")   = a1_v,
    Named("a2")   = a2_v
  );
}

// Read the FAM file, return IDs
// [[Rcpp::export]]
DataFrame read_fam_cpp(std::string fam_file) {
  std::ifstream f(fam_file);
  if (!f.is_open()) stop("Cannot open FAM file: " + fam_file);

  std::vector<std::string> fid, iid, pid, mid, sex_s;
  std::vector<double> pheno;
  std::string s_fid, s_iid, s_pid, s_mid, s_sex;
  double s_pheno;
  while (f >> s_fid >> s_iid >> s_pid >> s_mid >> s_sex >> s_pheno) {
    fid.push_back(s_fid);  iid.push_back(s_iid);
    pid.push_back(s_pid);  mid.push_back(s_mid);
    sex_s.push_back(s_sex); pheno.push_back(s_pheno);
  }
  f.close();
  return DataFrame::create(
    Named("FID") = fid, Named("IID") = iid,
    Named("PID") = pid, Named("MID") = mid,
    Named("sex") = sex_s, Named("pheno") = pheno
  );
}

// Read a batch of SNPs from a PLINK BED file
// snp_start: 0-indexed first SNP, snp_count: number of SNPs to read
// Returns n_indiv x snp_count matrix of dosages {0,1,2,NA}
// [[Rcpp::export]]
NumericMatrix read_bed_batch_cpp(std::string bed_file,
                                 int n_indiv,
                                 int snp_start,
                                 int snp_count) {
  std::ifstream f(bed_file, std::ios::binary);
  if (!f.is_open()) stop("Cannot open BED file: " + bed_file);

  // Check magic bytes and mode byte
  unsigned char magic[3];
  f.read(reinterpret_cast<char*>(magic), 3);
  if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
    f.close();
    stop("BED file does not appear to be SNP-major format.");
  }

  // Each SNP occupies ceil(n_indiv/4) bytes
  int bytes_per_snp = (n_indiv + 3) / 4;

  // Lookup: 2-bit code -> dosage {0,1,2,NA}
  // PLINK encoding: 00=hom A1, 01=missing, 10=het, 11=hom A2
  // Dosage convention: A2 = minor allele -> count of A2 alleles
  // 00 -> 0 (hom A1), 10 -> 1 (het), 11 -> 2 (hom A2), 01 -> NA
  static const double code_to_dose[4] = {0.0, NA_REAL, 1.0, 2.0};

  NumericMatrix G(n_indiv, snp_count);

  for (int j = 0; j < snp_count; j++) {
    // Seek to this SNP's data
    std::streampos offset = 3 + (std::streampos)(snp_start + j) * bytes_per_snp;
    f.seekg(offset);

    std::vector<unsigned char> buf(bytes_per_snp);
    f.read(reinterpret_cast<char*>(buf.data()), bytes_per_snp);

    for (int i = 0; i < n_indiv; i++) {
      int byte_idx = i / 4;
      int bit_idx  = (i % 4) * 2;
      int code     = (buf[byte_idx] >> bit_idx) & 0x03;
      G(i, j)      = code_to_dose[code];
    }
  }
  f.close();
  return G;
}

// Convenience: read single SNP
// [[Rcpp::export]]
NumericVector read_bed_snp_cpp(std::string bed_file, int n_indiv,
                               int snp_idx) {
  NumericMatrix G = read_bed_batch_cpp(bed_file, n_indiv, snp_idx, 1);
  return G(_, 0);
}

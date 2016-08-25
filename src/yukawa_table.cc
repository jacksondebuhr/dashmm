// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================

/// \file src/yukawa_table.cc
/// \brief Implementation of precomputed tables for Yukawa

#include "builtins/yukawa_table.h"

namespace dashmm {

const int YukawaTable::maxlev = 21; 
std::unique_ptr<YukawaTable> builtin_yukawa_table_; 

YukawaTable::YukawaTable(int n_digits, double size, double lambda) {
  const int maxlev = 21; 
  int p_table[] = {0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int s_table[] = {0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  lambda_ = lambda; 
  size_ = size; 
  scale_ = (lambda * size > 1.0 ? 1.0 / size : lambda) * size; 
  p_ = p_table[n_digits]; 
  s_ = s_table[n_digits]; 

  generate_sqf(); 
  generate_scaled_wigner_dmat(); 
  generate_m2m(); 
  generate_l2l(); 

  x_ = new double[s_]; 
  w_ = new double[s_]; 
  m_ = new int[s_ * (maxlev + 1)]; 
  sm_ = new int[(s_ + 1) * (maxlev + 1)]; 
  nexp_ = new int[maxlev + 1]; 
  f_ = new int[s_]; 
  smf_ = new int[(s_ + 1) * (maxlev + 1)]; 
  ealphaj_ = new dcomplex_t *[maxlev + 1]; 
  xs_ = new dcomplex_t *[maxlev + 1]; 
  ys_ = new dcomplex_t *[maxlev + 1]; 
  zs_ = new double *[maxlev + 1]; 

  int *m0 = new int[s_]; 

  switch (n_digits) {
  case 3: 
    x_[0] = 0.99273996739714473469540223504736787e-01;
    x_[1] = 0.47725674637049431137114652301534079e+00;
    x_[2] = 0.10553366138218296388373573790886439e+01;
    x_[3] = 0.17675934335400844688024335482623428e+01;
    x_[4] = 0.25734262935147067530294862081063911e+01;
    x_[5] = 0.34482433920158257478760788217186928e+01;
    x_[6] = 0.43768098355472631055818055756390095e+01;
    x_[7] = 0.53489575720546005399569367000367492e+01;
    x_[8] = 0.63576578531337464283978988532908261e+01; 

    w_[0] = 0.24776441819008371281185532097879332e+00;
    w_[1] = 0.49188566500464336872511239562300034e+00;
    w_[2] = 0.65378749137677805158830324216978624e+00; 
    w_[3] = 0.76433038408784093054038066838984378e+00;
    w_[4] = 0.84376180565628111640563702167128213e+00;
    w_[5] = 0.90445883985098263213586733400006779e+00;
    w_[6] = 0.95378613136833456653818075210438110e+00;
    w_[7] = 0.99670261613218547047665651916759089e+00;
    w_[8] = 0.10429422730252668749528766056755558e+01;
    
    f_[0] = 2;
    f_[1] = 4;
    f_[2] = 4; 
    f_[3] = 6;
    f_[4] = 6;
    f_[5] = 4;
    f_[6] = 6;
    f_[7] = 4;
    f_[8] = 2;  

    m0[0] = 4; 
    m0[1] = 8;
    m0[2] = 12;
    m0[3] = 16;
    m0[4] = 20; 
    m0[5] = 20; 
    m0[6] = 24; 
    m0[7] = 8; 
    m0[8] = 2;

    break; 
  default: 
    break;
  }

  for (int lev = 0; lev <= maxlev; ++lev) {
    double ld = lambda_ * size_ / pow(2, lev); 
    int *m = &m_[s_ * lev]; 
    int *sm = &sm_[(s_ + 1) * lev]; 
    int *smf = &smf_[(s_ + 1) * lev]; 

    for (int j = 0; j < s_; ++j) {
      double t1 = x_[j]; 
      double t2 = sqrt(t1 * t1 + 2 * t1 * ld); 
      int index = j;
      int max_mk = m0[j]; 
      for (int k = j + 1; k < s_; ++k) {
        if (x_[k] >= t2) {
          index = k; 
          break;
        } else {
          max_mk = (max_mk >= m0[k] ? max_mk : m0[k]);
        }
      }
      m[j] = (max_mk >= m0[index] ? max_mk : m0[index]);
    }

    sm[0] = 0; 
    smf[0] = 0; 
    for (int j = 1; j <= s_; ++j) {
      sm[j] = sm[j - 1] + m[j - 1] / 2; 
      smf[j] = smf[j - 1] + m[j - 1] * f_[j - 1] / 2;
    }
    nexp_[lev] = sm[s_]; 
  }

  for (int lev = 0; lev <= maxlev; ++lev) {
    int *smf = &smf_[(s_ + 1) * lev]; 
    int *m = &m_[s_ * lev]; 
    ealphaj_[lev] = new dcomplex_t[smf[s_]]; 
    
    int offset = 0; 
    for (int k = 0; k < s_; ++k) {
      double alpha = 2 * M_PI / m[k]; 
      for (int j = 1; j <= m[k] / 2; ++j) {
        double alpha_j = alpha * j; 
        for (int m = 1; m <= f_[k]; ++m) {
          ealphaj_[lev][offset++] = 
            dcomplex_t{cos(m * alpha_j), sin(m * alpha_j)}; 
        }
      }
    }
  } 

  for (int lev = 0; lev <= maxlev; ++lev) {
    xs_[lev] = new dcomplex_t[7 * nexp_[lev]]; 
    ys_[lev] = new dcomplex_t[7 * nexp_[lev]]; 
    zs_[lev] = new double[4 * nexp_[lev]]; 
    int *m = &m_[s_ * lev]; 
    double ld = lambda_ * size_ / pow(2, lev); 

    int offset = 0; 
    for (int k = 0; k < s_; ++k) {
      double alpha = 2 * M_PI / m[k];       
      for (int j = 1; j <= m[k] / 2; ++j) {
        double alphaj = j * alpha; 
        double arg = sqrt(x_[k] * x_[k] + 2 * x_[k] * ld) * cos(alphaj); 
        for (int d = -3; d <= 3; ++d) 
          xs_[lev][offset++] = dcomplex_t{cos(arg * d), sin(arg * d)};
      }
    }

    offset = 0; 
    for (int k = 0; k < s_; ++k) {
      double alpha = 2 * M_PI / m[k]; 
      for (int j = 1; j <= m[k] / 2; ++j) {
        double alphaj = j * alpha;
        double arg = sqrt(x_[k] * x_[k] + 2 * x_[k] * ld) * sin(alphaj);
        for (int d = -3; d <= 3; ++d) 
          ys_[lev][offset++] = dcomplex_t{cos(arg * d), sin(arg * d)};
      }
    }

    offset = 0; 
    for (int k = 0; k < s_; ++k) {
      for (int m = 0; m <= 3; ++m) {
        zs_[lev][offset++] = exp(-(x_[k] + ld) * m);
      }
    }
  }

  delete [] m0; 
}


YukawaTable::~YukawaTable() {
  delete [] sqf_; 
  for (auto it = dmat_plus_->begin(); it != dmat_plus_->end(); ++it) 
    delete [] it->second; 
  for (auto it = dmat_minus_->begin(); it != dmat_minus_->end(); ++it) 
    delete [] it->second; 
  delete dmat_plus_; 
  delete dmat_minus_; 
  delete [] m2m_; 
  delete [] l2l_;
  delete [] x_; 
  delete [] w_; 
  delete [] m_; 
  delete [] sm_; 
  delete [] nexp_; 
  delete [] f_; 
  delete [] smf_; 
  for (int i = 0; i <= maxlev; ++i) {
    delete [] ealphaj_[i]; 
    delete [] xs_[i]; 
    delete [] ys_[i];
    delete [] zs_[i];
  }
  delete [] ealphaj_; 
  delete [] xs_;
  delete [] ys_;
  delete [] zs_; 
}

void YukawaTable::generate_sqf() {
  sqf_ = new double[(p_ + 1) * (p_ + 2) / 2]; 
  double *temp = new double[2 * p_ + 1]; 

  temp[0] = 1.0; 
  for (int i = 1; i <= p_ * 2; ++i) 
    temp[i] = temp[i - 1] * i;

  int idx = 0; 
  for (int n = 0; n <= p_; ++n) {
    int n1 = 2 * n + 1; 
    for (int m = 0; m <= n; ++m) {
      sqf_[idx++] = n1 * temp[n - m] / temp[n + m];
    }
  }

  delete [] temp;
}

void YukawaTable::generate_scaled_wigner_dmat() {
  dmat_plus_ = new yukawa_map_t; 
  dmat_minus_ = new yukawa_map_t; 

  double cbeta[3] = {sqrt(3) / 3, -sqrt(3) / 3, 0}; 
  int nd = (p_ + 1) * (4 * p_ * p_ + 11 * p_ + 6) / 6; 
  for (int i = 0; i < 3; ++i) {
    double beta = acos(cbeta[i]); 
    double *dp_data = new double[nd]; 
    double *dm_data = new double[nd]; 
    generate_scaled_dmat_of_beta(beta, dp_data, dm_data); 
    (*dmat_plus_)[cbeta[i]] = dp_data;
    (*dmat_minus_)[cbeta[i]] = dm_data;
  }
}

void YukawaTable::generate_scaled_dmat_of_beta(double beta, double *dp, 
                                               double *dm) {
  double cbeta = cos(beta); 
  double sbeta = sin(beta); 
  double s2beta2 = (1 - cbeta) / 2; // sin^2(beta / 2)
  double c2beta2 = (1 + cbeta) / 2; // cos^2(beta / 2)

  // Set d_0^{0, 0} to 1 
  dp[0] = 1; 
  dm[0] = 1; 

  // Set d_1^{0, m} 
  dp[1] = -sbeta / sqrt(2); // d_1^{0, -1} 
  dp[2] = cbeta; // d_1^{0, 0}
  dp[3] = sbeta / sqrt(2); // d_1^{0, 1} 
  dm[1] = -dp[1]; 
  dm[2] = dp[2]; 
  dm[3] = -dp[3]; 

  // Set d_1^{1, m} 
  dp[4] = s2beta2; // d_1^{1, -1}
  dp[5] = -sbeta / sqrt(2); // d_1^{1, 0}
  dp[6] = c2beta2; // d_1^{1, 1}
  dm[4] = dp[4]; 
  dm[5] = -dp[5];
  dm[6] = dp[6];

  // Compute d_n^{0, m} for 2 <= n <= P
  for (int n = 2; n <= p_; ++n) {
    double *dpc = NULL, *dpp = NULL, *dmc = NULL; 
    int m; 
    // Get address of d_n^{0, 0} for angle beta, saved in dpc
    // Get address of d_{n - 1}^{0, 0} for angle beta, saved in dpp
    // Get address of d_n^{0, 0} for angle -beta, saved in dmc
    dpc = &dp[didx(n, 0, 0)]; 
    dpp = &dp[didx(n - 1, 0, 0)]; 
    dmc = &dm[didx(n, 0, 0)]; 

    // Compute d_n^{0, -n} 
    m = -n; 
    dpc[m] = -sbeta * 0.5 * sqrt((n - m) * (n - m - 1)) / n * dpp[m + 1]; 
    dmc[m] = dpc[m] * pow_m1(m); 

    // Compute d_n^{0, -n + 1}
    m = -n + 1; 
    dpc[m] = (-sbeta * 0.5 * sqrt((n - m) * (n - m - 1)) * dpp[m + 1] + 
              cbeta * sqrt((n - m) * (n + m)) * dpp[m]) / n; 
    dmc[m] = dpc[m] * pow_m1(m); 

    // Handle d_n^{0, m} 
    for (m = -n + 2; m <= n - 2; ++m) {
      dpc[m] = (-sbeta * 0.5 * sqrt((n - m) * (n - m - 1)) * dpp[m + 1] + 
                cbeta * sqrt((n - m) * (n + m)) * dpp[m] + 
                sbeta * 0.5 * sqrt((n + m) * (n + m - 1)) * dpp[m - 1]) / n; 
      dmc[m] = dpc[m] * pow_m1(m); 
    }

    // Handle d_n^{0, n - 1}
    m = n - 1; 
    dpc[m] = (cbeta * sqrt((n - m) * (n + m)) * dpp[m] +
              sbeta * 0.5 * sqrt((n + m) * (n + m - 1)) * dpp[m - 1]) / n;
    dmc[m] = dpc[m] * pow_m1(m); 

    // Handle d_n^{0, n} 
    m = n; 
    dpc[m] = sbeta * 0.5 * sqrt((n + m) * (n + m - 1)) / n * dpp[m - 1]; 
    dmc[m] = dpc[m] * pow_m1(m); 

    // Compute d_n^{mp, m} for 1 <= mp <= n 
    for (int mp = 1; mp <= n; ++mp) {
      // Get address of d_n^{mp, 0} for angle beta, saved in dpc
      // Get address of d_{n - 1}^{mp - 1, 0} for angle beta, saved in dpp
      // Get address of d_n^{mp, 0} for angle -beta, saved in dmc 
      dpc = &dp[didx(n, mp, 0)]; 
      dpp = &dp[didx(n - 1, mp - 1, 0)]; 
      dmc = &dm[didx(n, mp, 0)]; 

      double factor = 1.0 / sqrt((n + mp) * (n + mp - 1)); 

      // Compute d_n^{mp, -n}
      m = -n; 
      dpc[m] = s2beta2 * sqrt((n - m) * (n - m - 1)) * factor * dpp[m + 1]; 
      dmc[m] = dpc[m] * pow_m1(mp - m); 

      // Compute d_n^{mp, -n + 1}
      m = -n + 1; 
      dpc[m] = (s2beta2 * sqrt((n - m) * (n - m - 1)) * dpp[m + 1] -
                sbeta * sqrt((n + m) * (n - m)) * dpp[m]) * factor; 
      dmc[m] = dpc[m] * pow_m1(mp - m); 

      // Compute d_n^{mp, m}
      for (m = -n + 2; m <= n - 2; m++) {
        dpc[m] = (s2beta2 * sqrt((n - m) * (n - m - 1)) * dpp[m + 1] - 
                  sbeta * sqrt((n + m) * (n - m)) * dpp[m] + 
                  c2beta2 * sqrt((n + m) * (n + m - 1)) * dpp[m - 1]) * factor;
        dmc[m] = dpc[m] * pow_m1(mp - m); 
      }

      // Compute d_n^{mp, n - 1}
      m = n - 1; 
      dpc[m] = (-sbeta * sqrt((n + m) * (n - m)) * dpp[m] + 
                c2beta2 * sqrt((n + m) * (n + m - 1)) * dpp[m - 1]) * factor; 
      dmc[m] = dpc[m] * pow_m1(mp - m); 

      // Compute d_n^{mp, n}
      m = n; 
      dpc[m] = c2beta2 * sqrt((n + m) * (n + m - 1)) * factor * dpp[m - 1]; 
      dmc[m] = dpc[m] * pow_m1(mp - m); 
    }
  }

  // Scale 
  double factorial[2 * p_ + 1]; 
  factorial[0] = 1.0; 
  for (int i = 1; i <= 2 * p_; ++i) 
    factorial[i] = factorial[i - 1] * i; 

  int offset = 0; 
  for (int n = 0; n <= p_; ++n) {
    for (int mp = 0; mp <= n; ++mp) {
      for (int m = -n; m <= n; ++m) {
        double scale = sqrt(factorial[n - mp] * factorial[n + abs(m)] / 
                            factorial[n + mp] / factorial[n - abs(m)]); 
        dp[offset] *= scale; 
        dm[offset] *= scale;
        offset++;
      }
    }
  }
}

void YukawaTable::generate_m2m() {  
  m2m_ = new double[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * (maxlev + 1)]; 
  double *factorial = new double[2 * p_ + 1];
  double *bessel = new double[2 * p_ + 1]; 

  factorial[0] = 1.0; 
  for (int i = 1; i <= 2 * p_; ++i) 
    factorial[i] = factorial[i - 1] * i; 

  for (int lev = 0; lev <= maxlev; ++lev) {
    // lev refers to parent box
    
    // Compute shift distance rho
    double rho = sqrt(3) / 2 * size_ / pow(2, lev); 
    
    // Compute scaling factor 
    double sigma = scale_ / pow(2, lev); 

    // Compute scaled modified spherical bessel function
    bessel_in_scaled(2 * p_, lambda_ * rho, lambda_ * rho, bessel); 

    double *coeff = &m2m_[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * lev]; 

    for (int n = 0; n <= p_; ++n) {
      for (int m = 0; m <= n; ++m) {
        for (int np = m; np <= p_; ++np) {
          double temp = 0; 
          int bound = (n <= np ? n : np); 
          for (int k = m; k <= bound; ++k) {
            temp += pow(0.5, n + k) * pow(sigma, np - n)  * pow(-1, np + n) * 
              (2 * n + 1) * factorial[n - m] * factorial[np + m] * 
              factorial[2 * k] / factorial[k + m] / factorial[k - m] / 
              factorial[k] / factorial[n - k] / factorial[np - k] * 
              bessel[np + n - k] * pow(lambda_ * rho, np + n - 2 * k); 
          }
          coeff[sidx(n, m, np, p_)] = temp;
        }
      }
    }
  }

  delete [] factorial;
  delete [] bessel;
}

void YukawaTable::generate_l2l() {
  l2l_ = new double[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * (maxlev + 1)]; 
  double *factorial = new double[2 * p_ + 1]; 
  double *bessel = new double[2 * p_ + 1]; 

  factorial[0] = 1.0; 
  for (int i = 1; i <= 2 * p_; ++i) 
    factorial[i] = factorial[i - 1] * i; 

  for (int lev = 0; lev <= maxlev; ++lev) {
    // lev refers to the parent box 
    
    // Compute shift distance rho
    double rho = sqrt(3) / 4 * size_ / pow(2, lev); 

    // Compute scaling factor
    double sigma = scale_ / pow(2, lev); 

    // Compute scaled modified spherical bessel function
    bessel_in_scaled(2 * p_, lambda_ * rho, lambda_ * rho, bessel); 

    // Compute shifting coefficient
    double *coeff = &l2l_[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * lev]; 

    for (int n = 0; n <= p_; ++n) {
      for (int m = 0; m <= n; ++m) {
        for (int np = m; np <= p_; ++np) {
          double temp = 0; 
          int bound = (n <= np ? n : np); 
          for (int k = m; k <= bound; ++k) {
            temp += pow(0.5, n + k) * pow(sigma, n - np) * pow(-1, np + n) * 
              (2 * n + 1) * factorial[n - m] * factorial[np + m] * 
              factorial[2 * k] * bessel[n + np - k] *
              pow(lambda_ * rho, np + n - 2 * k) / factorial[k + m] / 
              factorial[k - m] / factorial[k] / factorial[n - k] / 
              factorial[np - k];
          }
          coeff[sidx(n, m, np, p_)] = temp;
        }
      }
    }
  }

  delete [] factorial; 
  delete [] bessel;
}

void update_yukawa_table(int n_digits, double size, double lambda) {
  if (builtin_yukawa_table_ == nullptr) 
    builtin_yukawa_table_ = 
      std::unique_ptr<YukawaTable>{new YukawaTable{n_digits, size, lambda}}; 
}

} // namespace dashmm

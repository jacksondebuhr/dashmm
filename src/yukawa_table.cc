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
const double M_SQRTPI = 0.9189385332046727417803297L;
uYukawaTable builtin_yukawa_table_; 

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

void legendre_Pnm(int nb, double x, double *P) {
  double u = -sqrt(1.0 - x * x); 
  P[midx(0, 0)] = 1.0; 
  for (int n = 1; n <= nb; ++n) 
    P[midx(n, n)] = P[midx(n - 1, n - 1)] * u * (2 * n - 1); 

  for (int n = 0; n < nb; ++n) 
    P[midx(n + 1, n)] = P[midx(n, n)] * x * (2 * n + 1); 

  for (int m = 0; m <= nb; ++m) {
    for (int n = m + 2; n <= nb; ++n) {
      P[midx(n, m)] = ((2.0 * n - 1) * x * P[midx(n - 1, m)] - 
                         (n + m - 1) * P[midx(n - 2, m)]) / (n - m);
    }
  }
}

void legendre_Pnm_gt1_scaled(int nb, double x, double scale, double *P) {
  double v = scale * x; 
  double w = scale * scale; 
  double u = sqrt(x * x - 1.0) * scale; 

  P[midx(0, 0)] = 1.0; 
  for (int n = 1; n <= nb; ++n) 
    P[midx(n, n)] = P[midx(n - 1, n - 1)] * u * (2 * n - 1); 

  for (int n = 0; n < nb; ++n) 
    P[midx(n + 1, n)] = P[midx(n, n)] * (2 * n + 1) * v; 

  for (int m = 0; m <= nb; ++m) {
    for (int n = m + 2; n <= nb; ++n) {
      P[midx(n, m)] = ((2.0 * n - 1) * v * P[midx(n - 1, m)] 
                       - (n + m - 1) * w * P[midx(n - 2, m)]) / (n - m); 
    }
  }
}

double Gamma(double x) {
  // the largest argument for which gamma(x) is representable
  const double xbig = 171.624; 

  // the smallest positive floating-point number such that 1/xminin is
  // machine representable  
  const double xminin = 2.23e-308;

  // the smallest positive floating-point number such that 1 + eps > 1
  const double eps = 2.22e-16;

  // the largest machine representable floating-point number 
  const double xinf = 1.79e308; 

  // numerator and denominator coefficients for rational minimax
  // approximation over (1, 2) 
  const double p[8] = {-1.71618513886549492533811e+0,
                       2.47656508055759199108314e+1,
                       -3.79804256470945635097577e+2,
                       6.29331155312818442661052e+2,
                       8.66966202790413211295064e+2,
                       -3.14512729688483675254357e+4,
                       -3.61444134186911729807069e+4,
                       6.64561438202405440627855e+4}; 

  const double q[8] = {-3.08402300119738975254353e+1,
                       3.15350626979604161529144e+2,
                       -1.01515636749021914166146e+3,
                       -3.10777167157231109440444e+3,
                       2.25381184209801510330112e+4,
                       4.75584627752788110767815e+3,
                       -1.34659959864969306392456e+5,
                       -1.15132259675553483497211e+5}; 

  // coefficients for minimax approximation over (12, inf) 
  const double c[7] = {-1.910444077728e-03,
                       8.4171387781295e-04,
                       -5.952379913043012e-04,
                       7.93650793500350248e-04,
                       -2.777777777777681622553e-03,
                       8.333333333333333331554247e-02,
                       5.7083835261e-03}; 

  bool parity = false; 
  double y, fact, res; 

  fact = 1; 
  y = x; 

  if (y <= 0.0) {
    // argument is negative
    y = -x; 
    double y1 = floor(y); 
    res = y - y1; 
    if (res != 0.0) {
      if (y1 != floor(y1 * 0.5) * 2) 
        parity = true;
      fact = -M_PI / sin(M_PI * res);
      y += 1.0;
    } else {
      return xinf; 
    }
  }

  // argument is positive
  if (y < eps) {
    if (y >= xminin) {
      res = 1.0 / y;
    } else {
      return xinf;
    } 
  } else if (y < 12.0) {
    double y1 = y; 
    double z; 
    int n; 
    if (y < 1.0) {
      z = y; 
      y += 1.0;
    } else {
      n = floor(y) - 1; 
      y -= n; 
      z = y - 1.0; 
    } 

    double xnum = 0.0; 
    double xden = 1.0; 
    
    for (int i = 0; i < 8; ++i) {
      xnum = (xnum + p[i]) * z; 
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0; 
    
    if (y1 < y) {
      res = res / y; 
    } else if (y1 > y) {
      for (int i = 0; i < n; ++i) {
        res *= y; 
        y += 1.0;
      }
    }
  } else {
    if (y <= xbig) {
      double ysq = y * y; 
      double sum = c[6]; 
      for (int i = 0; i < 6; ++i) {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + M_SQRTPI;
      sum += (y - 0.5) * log(y); 
      res = exp(sum); 
    } else {
      return xinf;
    }
  }

  if (parity) 
    res = -res; 
  if (fact != 1.0) 
    res = fact / res; 
  return res; 
}

int bessel_In(int nb, double alpha, double x, int ize, double *B) {
  // decimal significance desired
  const int nsig = 16; 
  
  // upper limit on the magnitude of x when ize is 1
  const double exparg = 709.0; 
  
  // upper limit on the magnitude of x when ize is 2
  const double xlarge = 1e4; 

  // 10^k, where k is the largest integer such that enten is
  // machine-representable in working precision
  const double enten = 1e308; 

  // 10^nsig
  const double ensig = 1e16; 

  // 10^(-k) for the smallest integer k such that k >= nsig / 4
  const double rtnsig = 1e-4; 

  // smallest abs(x) such that x / 4 does not underflow
  const double enmten = 8.9e-308; 

  // Return if nb, alpha, x, or ize are out of range
  if (nb <= 0 || x < 0 || alpha < 0 || alpha >= 1 ||
      (ize == 1 && x > exparg) || 
      (ize == 2 && x > xlarge)) {
    return -1;
  }

  int ncalc = nb; 
  double tempa, tempb, tempc; 
  double empal, emp2al, halfx, tover, test;
  double p, plast, pold, psave, psavel, sum, en, em; 
  int magx, nbmx, n, nstart, nend; 

  if (x >= rtnsig) {
    magx = (int) x; 
    nbmx = nb - magx; 
    n = magx + 1; 
    en = (double)(n + n) + (alpha + alpha); 
    plast = 1.0; 
    p = en / x; 

    // Calculate general significance test
    test = ensig + ensig; 
    test = (2 * magx > 5 * nsig ? sqrt(test * p) : 
            test / pow(1.585, magx)); 

    bool stepin = true; 

    if (nbmx >= 3) {
      // Calculate p-sequence until n = nb - 1. Check for possible overflow
      tover = enten / ensig; 
      nstart = magx + 2; 
      nend = nb - 1; 
      for (int k = nstart; k <= nend; ++k) {
        n = k; 
        en += 2.0; 
        pold = plast; 
        plast = p; 
        p = en * plast / k + pold; 
        if (p > tover) {
          // To avoid overflow, divide p-sequence by tover. Calculate p-sequence
          // until abs(p) > 1 
          tover = enten; 
          p /= tover; 
          plast /= tover;
          psave = p; 
          psavel = plast; 
          nstart = n + 1; 
          do {
            n += 1; 
            en += 2.0; 
            pold = plast; 
            plast = p; 
            p = en * plast / x + pold; 
          } while (p <= 1.0); 
          tempb = en / x; 

          // Calculate backward test, and find ncalc, the highest n such that
          // the test is passed 
          test = pold * plast / ensig; 
          test = test * (0.5 - 0.5 / (tempb * tempb)); 
          p = plast * tover; 
          n -= 1; 
          en -= 2.0; 
          nend = (nb <= n ? nb : n); 
          for (int ell = nstart; ell <= nend; ++ell) {
            ncalc = ell; 
            pold = psavel; 
            psavel = psave; 
            psave = en * psavel / x + pold; 
            if (psave * psavel > test) 
              break; 
          }
          if (ncalc != nend) 
            ncalc--; 

          stepin = false; 
          break;
        }      
      }
      
      if (stepin) { 
        n = nend; 
        en = (double)(n + n) + (alpha + alpha);       
        // Calculate special significance test for nbmx >= 3
        test = fmax(test, sqrt(plast * ensig) * sqrt(p + p)); 
      }
    }

    if (stepin) {
      // Calculate p-sequence 
      do {
        n += 1; 
        en += 2.0; 
        pold = plast; 
        plast = p; 
        p = en * plast / x + pold; 
      } while (p < test); 
    }

    // Initialize the backward recursion and the normalization sum 
    n += 1; 
    en += 2.0; 
    tempb = 0.0; 
    tempa = 1.0 / p; 
    em = (double) n - 1.0; 
    empal = em + alpha; 
    emp2al = (em - 1.0) + (alpha + alpha); 
    sum = tempa * empal * emp2al / em; 
    nend = n - nb; 

    if (nend < 0) {
      // Store B[n - 1] and set higher orders to 0.0
      B[n - 1] = tempa; 
      nend = -nend; 
      for (int ell = 0; ell < nend; ++ell) 
        B[n + ell] = 0.0;
    } else {
      // Recur backward via difference equation, calculating (but not storing
      // b[n - 1]), until n = nb
      for (int ell = 0; ell < nend; ++ell) {
        n -= 1; 
        en -= 2.0; 
        tempc = tempb; 
        tempb = tempa; 
        tempa = en * tempb / x + tempc; 
        em -= 1.0; 
        emp2al -= 1.0; 
        if (n == 1) 
          break; 
        if (n == 2) 
          emp2al = 1.0; 
        empal -= 1.0; 
        sum = (sum + tempa * empal) * emp2al / em; 
      }

      // Store B[nb - 1] 
      B[n - 1] = tempa; 
      if (nb <= 1) {
        sum = (sum + sum) + tempa; 
        if (alpha != 0.0)  // statement 230
          sum = sum * Gamma(alpha + 1.0) * pow(x * 0.5, -alpha); 
        if (ize == 1) 
          sum *= exp(-x); 
        tempa = enmten; 
        if (sum > 1.0) 
          tempa *= sum; 
        
        for (int n = 0; n < nb; ++n) {
          if (B[n] < tempa) 
            B[n] = 0.0; 
          B[n] /= sum;
        }    

        return ncalc; 
      }

      // Calculate and store B[nb - 2] 
      n -= 1; 
      en -= 2.0; 
      B[n - 1] = (en * tempa) / x + tempb; 
      if (n == 1) {
        sum = (sum + sum) + B[0]; 
        if (alpha != 0.0)  // statement 230
          sum = sum * Gamma(alpha + 1.0) * pow(x * 0.5, -alpha); 
        if (ize == 1) 
          sum *= exp(-x); 
        tempa = enmten; 
        if (sum > 1.0) 
          tempa *= sum; 
        
        for (int n = 0; n < nb; ++n) {
          if (B[n] < tempa) 
            B[n] = 0.0; 
          B[n] /= sum;
        }    
        
        return ncalc;        
      } 
      em -= 1.0; 
      emp2al -= 1.0; 
      if (n == 2) 
        emp2al = 1.0; 
      empal -= 1.0; 
      sum = (sum + B[n - 1] * empal) * emp2al / em;          
    }
    
    nend = n - 2; 

    if (nend > 0) {
      // Calculate via difference equation and store B[n - 1], until n = 2
      for (int ell = 0; ell < nend; ++ell) {
        n -= 1; 
        en -= 2.0; 
        B[n - 1] = (en * B[n]) / x + B[n + 1]; 
        em -= 1.0; 
        emp2al -= 1.0; 
        if (n == 2) 
          emp2al = 1.0; 
        empal -= 1.0; 
        sum = (sum + B[n - 1] * empal) * emp2al / em;
      }
    }

    // Calculate B[0]
    B[0] = 2.0 * empal * B[1] / x + B[2]; 
    sum = (sum + sum) + B[0];  

    // Normalize, divide all B[n] by sum 
    if (alpha != 0.0)  // statement 230
      sum = sum * Gamma(alpha + 1.0) * pow(x * 0.5, -alpha); 
    if (ize == 1) 
      sum *= exp(-x); 
    tempa = enmten; 
    if (sum > 1.0) 
      tempa *= sum; 

    for (int n = 0; n < nb; ++n) {
      if (B[n] < tempa) 
        B[n] = 0.0; 
      B[n] /= sum;
    }    
  } else {
    empal = 1.0 + alpha; 
    halfx = (x > enmten ? x * 0.5 : 0.0); 
    tempa = (alpha != 0 ? pow(halfx, alpha) / Gamma(empal) : 1.0); 
    if (ize == 2) 
      tempa *= exp(-x); 
    tempb = ((x + 1.0) > 1.0 ? halfx * halfx : 0.0); 
    B[0] = tempa + tempa * tempb / empal; 
    if (x != 0 && B[0] == 0.0) 
      ncalc = 0; 

    if (nb > 1) {
      if (x == 0.0) {
        for (int i = 1; i < nb; ++i) 
          B[i] = 0.0;
      } else {
        // Calculate higher order functions
        tempc = halfx; 
        tover = (tempb != 0 ? enmten / tempb : 
                 (enmten + enmten) / x); 

        for (int n = 1; n < nb; ++n) {
          tempa /= empal; 
          empal += 1.0; 
          tempa *= tempc; 
          if (tempa <= tover * empal) 
            tempa = 0.0; 
          B[n] = tempa + tempa * tempb / empal; 
          if (B[n] == 0 && ncalc > n)
            ncalc = n - 1; 
        }
      }
    }
  }  

  return ncalc;     
}

void bessel_in_scaled(int nb, double x, double scale, double *B) {
  const double ensig = 1e-4; 
  const double enmten = 1e-300; 
  
  if (x <= ensig) {
    // Use 2-term Taylor expansion
    double scale_x = x / scale; 
    double t1 = 1.0; 
    double t2 = 0.5 * x * x; 
    B[0] = t1 * (1.0 + t2 / 3.0); 
    for (int i = 1; i <= nb; ++i) {
      t1 = t1 * scale_x / (2 * i + 1); 
      if (t1 <= enmten) 
        t1 = 0.0; 
      B[i] = t1 * (1.0 + t2 / (2 * i + 3));
    }
  } else if (x > 1.0e2) {
    for (int  i = 0; i <= nb; ++i) 
      B[i] = 0.0;
  } else {
    double factor = sqrt(M_PI_2 / x); 
    
    assert(bessel_In(nb + 1, 0.5, x, 1, B) == (nb + 1)); 
    
    for (int i = 0; i <= nb; ++i) {
      B[i] *= factor; 
      factor /= scale; 
      if (fabs(B[i]) <= enmten) 
        factor = 0.0; 
    }
  }
}

void bessel_kn_scaled(int nb, double x, double scale, double *B) {
  const double xmin = 4.46e-308;
  const double xinf = 1.79e308;
  const double xlarge = 1e8; 
  int ncalc = 0; 
  
  if (nb >= 0 && x >= xmin && x < xlarge) {
    double ex = x; 
    double p = exp(-ex) / ex * M_PI_2; 
    B[0] = p; 
    B[1] = p * scale * (1.0 + 1.0 / ex); 
    ncalc = 1; 
    double u1 = scale / ex; 
    double u2 = scale * scale; 
    for (int n = 2; n <= nb; ++n) {
      if (fabs(B[n - 1]) * u1 >= xinf / (2 * n - 1)) 
        break;
      B[n] = (2 * n - 1) * u1 * B[n - 1] + u2 * B[n - 2]; 
      ncalc++; 
    }

    for (int n = ncalc + 1; n <= nb; ++n) 
      B[n] = 0.0; 
  } else {
    B[0] = 0.0; 
    ncalc = (nb + 1 <= 0 ? nb + 1 : 0) - 1; 
  }
         
  assert(ncalc == nb); 
}

} // namespace dashmm

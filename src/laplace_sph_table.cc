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


/// \file src/laplace_sph_table.cc
/// \brief Implementation of precomputed tables for LaplaceSPH


#include "builtins/laplace_sph_table.h"


namespace dashmm {


std::map<int, uLaplaceSPHTable> builtin_laplace_table_;


void legendre_Plm(int n, double x, double *P) {
  double u = -sqrt(1.0 - x * x);
  P[midx(0, 0)] = 1.0;
  for (int i = 1; i <= n; i++)
    P[midx(i, i)] = P[midx(i - 1, i - 1)] * u * (2 * i - 1);

  for (int i = 0; i < n; i++)
    P[midx(i + 1, i)] = P[midx(i, i)] * x * (2 * i + 1);

  for (int m = 0; m <= n; m++) {
    for (int ell = m + 2; ell <= n; ell++) {
      P[midx(ell, m)] = ((2.0 * ell - 1) * x * P[midx(ell - 1, m)] -
                         (ell + m - 1) * P[midx(ell - 2, m)]) / (ell - m);
    }
  }
}


LaplaceSPHTable::LaplaceSPHTable(int n_digits) {
  int expan_length[] = {0, 4, 7, 9, 13, 16, 18, 23, 26, 29,
                        33, 36, 40, 43, 46};
  p_ = expan_length[n_digits];
  sqf_ = generate_sqf();
  sqbinom_ = generate_sqbinom();
  dmat_plus_ = new laplace_map_t;
  dmat_minus_ = new laplace_map_t;
  generate_wigner_dmatrix(dmat_plus_, dmat_minus_);
}


LaplaceSPHTable::~LaplaceSPHTable() {
  delete [] sqf_;
  delete [] sqbinom_;
  for (auto it = dmat_plus_->begin(); it != dmat_plus_->end(); ++it)
    delete [] it->second;
  for (auto it = dmat_minus_->begin(); it != dmat_minus_->end(); ++it)
    delete [] it->second;
  delete dmat_plus_;
  delete dmat_minus_;
}


double *LaplaceSPHTable::generate_sqf() {
  double *sqf = new double[2 * p_ + 1];
  sqf[0] = 1.0;
  for (int i = 1; i <= p_ * 2; ++i)
    sqf[i] = sqf[i - 1] * sqrt(i);
  return sqf;
}


double *LaplaceSPHTable::generate_sqbinom() {
  int N = 2 * p_;
  int total = (N + 1) * (N + 2) / 2;
  double *sqbinom = new double[total];

  // Set the first three binomial coefficients
  sqbinom[0] = 1; // binom(0, 0)
  sqbinom[1] = 1; // binom(1, 0)
  sqbinom[2] = 1; // binom(1, 1)

  // Compute the rest binomial coefficients
  for (int n = 2; n <= N; ++n) {
    // Get address of binom(n, 0) in current
    // Get address of binom(n - 1, 0) in previous
    double *current = &sqbinom[midx(n, 0)];
    double *previous = &sqbinom[midx(n - 1, 0)];

    current[0] = 1; // binom(n, 0);
    for (int m = 1; m < n; ++m) {
      current[m] = previous[m] + previous[m - 1];
    }
    current[n] = 1; // binom(n, n)
  }

  // Compute the square root of the binomial coefficients
  for (int i = 0; i < total; ++i) {
    sqbinom[i] = sqrt(sqbinom[i]);
  }

  return sqbinom;
}


void LaplaceSPHTable::generate_wigner_dmatrix(laplace_map_t *&dp,
                                              laplace_map_t *&dm) {
  double cbeta[24] = {1.0 / sqrt(5.0), 1.0 / sqrt(6.0), 1.0 / sqrt(9.0),
                      1.0 / sqrt(10.0), 1.0 / sqrt(11.0), 1.0 / sqrt(14.0),
                      1.0 / sqrt(19.0), 2.0 / sqrt(5.0), sqrt(2.0 / 3.0),
                      sqrt(1.0 / 2.0), sqrt(4.0 / 9.0), sqrt(1.0 / 3.0),
                      sqrt(4.0 / 13.0), sqrt(2.0 / 7.0), sqrt(4.0 / 17.0),
                      sqrt(2.0 / 11.0), sqrt(9.0 / 10.0), sqrt(9.0 / 11.0),
                      sqrt(9.0 / 13.0), sqrt(9.0 / 14.0), sqrt(9.0 / 17.0),
                      sqrt(9.0 / 19.0), sqrt(9.0 / 22.0), 0.0};
  int nd = (p_ + 1) * (4 * p_ * p_ + 11 * p_ + 6) / 6;
  for (int i = 0; i < 24; ++i) {
    double beta = acos(cbeta[i]);
    double *dp_data = new double[nd];
    double *dm_data = new double[nd];
    generate_dmatrix_of_beta(beta, dp_data, dm_data);
    (*dp)[cbeta[i]] = dp_data;
    (*dm)[cbeta[i]] = dm_data;
  }

  for (int i = 0; i < 23; i++) {
    double beta = acos(-cbeta[i]);
    double *dp_data = new double[nd];
    double *dm_data = new double[nd];
    generate_dmatrix_of_beta(beta, dp_data, dm_data);
    (*dp)[-cbeta[i]] = dp_data;
    (*dm)[-cbeta[i]] = dm_data;
  }
}


void LaplaceSPHTable::generate_dmatrix_of_beta(double beta,
                                               double *dp, double *dm) {
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
    for (int mp = 1; mp <= n; mp++) {
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
      for (m = -n + 2; m <= n - 2; ++m) {
        dpc[m] = (s2beta2 * sqrt((n - m) * (n - m - 1)) * dpp[m + 1]
                  + c2beta2 * sqrt((n + m) * (n + m - 1)) * dpp[m - 1]
                  - sbeta * sqrt((n + m) * (n - m)) * dpp[m]) * factor;
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
}


} // namespace dashmm

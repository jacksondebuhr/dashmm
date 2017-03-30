// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================

/// \file src/helmholtz_table.cc
/// \brief Implementation of precomputed tables for Helmholtz

#include "builtins/helmholtz_table.h"

namespace dashmm {

const int HelmholtzTable::maxlev = 21;
std::unique_ptr<HelmholtzTable> builtin_helmholtz_table_;

HelmholtzTable::HelmholtzTable(int n_digits, double size, double omega) {
  int p_table[] = {0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int evan_table[] = {0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int prop_table[] = {0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  omega_ = omega;
  size_ = size;
  scale_ = (omega * size > 1.0 ? 1.0  / size : omega) * size;
  p_ = p_table[n_digits];
  s_e_ = evan_table[n_digits];
  s_p_ = prop_table[n_digits];

  generate_sqf();
  generate_scaled_wigner_dmat();
  generate_m2m();
  generate_l2l();

  x_e_ = new double[s_e_]();
  w_e_ = new double[s_e_]();
  m_e_ = new int[s_e_ * (maxlev + 1)]();
  sm_e_ = new int[(s_e_ + 1) * (maxlev + 1)]();
  n_e_ = new int[maxlev + 1]();
  f_e_ = new int[s_e_]();
  smf_e_ = new int[(s_e_ + 1) * (maxlev + 1)]();
  ealphaj_e_ = new dcomplex_t *[maxlev + 1];
  xs_e_ = new dcomplex_t *[maxlev + 1];
  ys_e_ = new dcomplex_t *[maxlev + 1];
  zs_e_ = new double *[maxlev + 1];

  x_p_ = new double[s_p_]();
  w_p_ = new double[s_p_]();
  m_p_ = new int[s_p_ * (maxlev + 1)]();
  sm_p_ = new int[(s_p_ + 1) * (maxlev + 1)]();
  n_p_ = new int[maxlev + 1]();
  f_p_ = new int[s_p_]();
  smf_p_ = new int[(s_p_ + 1) * (maxlev + 1)]();
  ealphaj_p_ = new dcomplex_t *[maxlev + 1];
  xs_p_ = new dcomplex_t *[maxlev + 1];
  ys_p_ = new dcomplex_t *[maxlev + 1];
  zs_p_ = new dcomplex_t *[maxlev + 1];

  int *m0e = new int[s_e_];
  int *m0p = new int[s_p_];
  gaussq(s_p_);

  switch (n_digits) {
  case 3:
    x_e_[0] = 0.99273996739714473469540223504736787e-01;
    x_e_[1] = 0.47725674637049431137114652301534079e+00;
    x_e_[2] = 0.10553366138218296388373573790886439e+01;
    x_e_[3] = 0.17675934335400844688024335482623428e+01;
    x_e_[4] = 0.25734262935147067530294862081063911e+01;
    x_e_[5] = 0.34482433920158257478760788217186928e+01;
    x_e_[6] = 0.43768098355472631055818055756390095e+01;
    x_e_[7] = 0.53489575720546005399569367000367492e+01;
    x_e_[8] = 0.63576578531337464283978988532908261e+01;

    w_e_[0] = 0.24776441819008371281185532097879332e+00;
    w_e_[1] = 0.49188566500464336872511239562300034e+00;
    w_e_[2] = 0.65378749137677805158830324216978624e+00;
    w_e_[3] = 0.76433038408784093054038066838984378e+00;
    w_e_[4] = 0.84376180565628111640563702167128213e+00;
    w_e_[5] = 0.90445883985098263213586733400006779e+00;
    w_e_[6] = 0.95378613136833456653818075210438110e+00;
    w_e_[7] = 0.99670261613218547047665651916759089e+00;
    w_e_[8] = 0.10429422730252668749528766056755558e+01;

    f_e_[0] = 2;
    f_e_[1] = 3;
    f_e_[2] = 4;
    f_e_[3] = 5;
    f_e_[4] = 5;
    f_e_[5] = 4;
    f_e_[6] = 5;
    f_e_[7] = 4;
    f_e_[8] = 1;

    m0e[0] = 4;
    m0e[1] = 8;
    m0e[2] = 12;
    m0e[3] = 16;
    m0e[4] = 20;
    m0e[5] = 20;
    m0e[6] = 24;
    m0e[7] = 8;
    m0e[8] = 2;

    f_p_[0] = 10;
    f_p_[1] = 10;
    f_p_[2] = 10;
    f_p_[3] = 10;
    f_p_[4] = 10;
    f_p_[5] = 10;
    f_p_[6] = 10;
    f_p_[7] = 10;
    f_p_[8] = 10;
    f_p_[9] = 10;

    m0p[0] = 16;
    m0p[1] = 16;
    m0p[2] = 16;
    m0p[3] = 16;
    m0p[4] = 16;
    m0p[5] = 16;
    m0p[6] = 16;
    m0p[7] = 16;
    m0p[8] = 16;
    m0p[9] = 16;

    break;
  default:
    break;
  }

  for (int lev = 0; lev <= maxlev; ++lev) {
    double wd = omega_ * size_ / pow(2, lev);
    int C1, C2;
    if (wd <= 0.125) {
      C1 = 0;
      C2 = 0;
    } else if (omega <= 0.25) {
      C1 = 6;
      C2 = 4;
    } else if (omega <= 6) {
      C1 = 20;
      C2 = 8;
    } else if (omega <= 8) {
      C1 = 18;
      C2 = 10;
    } else if (omega <= 10) {
      C1 = 22;
      C2 = 12;
    } else {
      fprintf(stderr, "Exceeds low-frequency regime\n"); 
      exit(-1); 
    }

    int *m_e = &m_e_[s_e_ * lev];
    int *sm_e = &sm_e_[(s_e_ + 1) * lev];
    int *smf_e = &smf_e_[(s_e_ + 1) * lev];
    int *m_p = &m_p_[s_p_ * lev];
    int *sm_p = &sm_p_[(s_p_ + 1) * lev];
    int *smf_p = &smf_p_[(s_p_ + 1) * lev];

    if (wd > 1e-4) {
      for (int j = 0; j < s_e_; ++j) {
        double t1 = x_e_[j]; 
        double t2 = sqrt(t1 * t1 + wd * wd); 
        int index = j; 
        int max_mk = m0e[j]; 
        for (int k = j + 1; k < s_e_; ++k) {
          if (x_e_[k] >= t2) {
            index = k;
            break;
          } else {
            max_mk = (max_mk >= m0e[k] ? max_mk : m0e[k]);
          }
        }
        m_e[j] = (max_mk >= m0e[index] ? max_mk : m0e[index]) + C1;
      }
    } else {
      for (int j = 0; j < s_e_; ++j)
        m_e[j] = m0e[j];
    }

    for (int j = 0; j < s_p_; ++j)
      m_p[j] = m0p[j] + C2;

    sm_e[0] = 0;
    smf_e[0] = 0;
    for (int j = 1; j <= s_e_; ++j) {
      sm_e[j] = sm_e[j - 1] + m_e[j - 1] / 2;
      smf_e[j] = smf_e[j - 1] + m_e[j - 1] * f_e_[j - 1] / 2;
    }
    n_e_[lev] = sm_e[s_e_];

    sm_p[0] = 0; 
    smf_p[0] = 0; 
    for (int j = 1; j <= s_p_; ++j) {
      sm_p[j] = sm_p[j - 1] + m_p[j - 1]; 
      smf_p[j] = smf_p[j - 1] + m_p[j - 1] * f_p_[j - 1]; 
    }

    n_p_[lev] = sm_p[s_p_];
  }

  for (int lev = 0; lev <= maxlev; ++lev) {
    int *smf_e = &smf_e_[(s_e_ + 1) * lev];
    int *m_e = &m_e_[s_e_ * lev];
    ealphaj_e_[lev] = new dcomplex_t[smf_e[s_e_]];

    int offset = 0;
    for (int k = 0; k < s_e_; ++k) {
      double alpha_k = 2 * M_PI / m_e[k];
      for (int j = 1; j <= m_e[k] / 2; ++j) {
        double alpha_kj = alpha_k * j;
        for (int m = 1; m <= f_e_[k]; ++m) {
          ealphaj_e_[lev][offset++] =
            dcomplex_t{cos(m * alpha_kj), sin(m * alpha_kj)};
        }
      }
    }

    int *smf_p = &smf_p_[(s_p_ + 1) * lev];
    int *m_p = &m_p_[s_p_ * lev];
    ealphaj_p_[lev] = new dcomplex_t[smf_p[s_p_]];

    offset = 0; 
    for (int k = 0; k < s_p_; ++k) {
      double alpha_k = 2 * M_PI / m_p[k]; 
      for (int j = 1; j <= m_p[k]; ++j) {
        double alpha_kj = alpha_k * j; 
        for (int m = 1; m <= f_p_[k]; ++m) {
          ealphaj_p_[lev][offset++] = 
            dcomplex_t{cos(m * alpha_kj), sin(m * alpha_kj)};
        }
      }
    }
  }

  for (int lev = 0; lev <= maxlev; ++lev) {
    xs_e_[lev] = new dcomplex_t[7 * n_e_[lev]];
    ys_e_[lev] = new dcomplex_t[7 * n_e_[lev]];
    zs_e_[lev] = new double[4 * s_e_];
    int *m_e = &m_e_[s_e_ * lev];
    double wd = omega_ * size_ / pow(2, lev);

    // Evanescent wave, x-direction shift coefficients
    int offset = 0;
    for (int k = 0; k < s_e_; ++k) {
      double alpha_k = 2 * M_PI / m_e[k];
      for (int j = 1; j <= m_e[k] / 2; ++j) {
        double alpha_kj = alpha_k * j;
        double arg = sqrt(x_e_[k] * x_e_[k] + wd * wd) * cos(alpha_kj);
        for (int d = -3; d <= 3; ++d) {
          xs_e_[lev][offset++] = dcomplex_t{cos(arg * d), sin(arg * d)};
        }
      }
    }

    // Evanescent wave, y-direction shift coefficients
    offset = 0;
    for (int k = 0; k < s_e_; ++k) {
      double alpha_k = 2 * M_PI / m_e[k];
      for (int j = 1; j <= m_e[k] / 2; ++j) {
        double alpha_kj = alpha_k * j;
        double arg = sqrt(x_e_[k] * x_e_[k] + wd * wd) * sin(alpha_kj);
        for (int d = -3; d <= 3; ++d) {
          ys_e_[lev][offset++] = dcomplex_t{cos(arg * d), sin(arg * d)};
        }
      }
    }

    // Evanescent wave, z-direction shift coefficients
    offset = 0;
    for (int k = 0; k < s_e_; ++k) {
      for (int m = 0; m <= 3; ++m) {
        zs_e_[lev][offset++] = exp(-x_e_[k] * m);
      }
    }

    xs_p_[lev] = new dcomplex_t[7 * n_p_[lev]];
    ys_p_[lev] = new dcomplex_t[7 * n_p_[lev]];
    zs_p_[lev] = new dcomplex_t[7 * s_p_];
    int *m_p = &m_p_[s_p_ * lev];

    // Propagating wave, x-direction shift coefficients
    offset = 0;
    for (int k = 0; k < s_p_; ++k) {
      double alpha_k = 2 * M_PI / m_p[k];
      for (int j = 1; j <= m_p[k]; ++j) {
        double alpha_kj = j * alpha_k; 
        double arg = wd * sin(x_p_[k]) * cos(alpha_kj); 
        for (int d = -3; d <= 3; ++d) {
          xs_p_[lev][offset++] = dcomplex_t{cos(arg * d), sin(arg * d)};
        }
      }
    }

    // Propagating wave, y-direction shift coefficients
    offset = 0;
    for (int k = 0; k < s_p_; ++k) {
      double alpha_k = 2 * M_PI / m_p[k];
      for (int j = 1; j <= m_p[k]; ++j) {
        double alpha_kj = j * alpha_k; 
        double arg = wd * sin(x_p_[k]) * sin(alpha_kj); 
        for (int d = -3; d <= 3; ++d) {
          ys_p_[lev][offset++] = dcomplex_t{cos(arg * d), sin(arg * d)};
        }
      }
    }

    // Propagating wave, z-direction shift coefficients
    offset = 0;
    for (int k = 0; k < s_p_; ++k) {
      double arg = wd * cos(x_p_[k]);
      for (int m = -3; m <= 3; ++m) {
        zs_p_[lev][offset++] = dcomplex_t{cos(arg * m), sin(arg * m)};
      }
    }
  }

  delete [] m0e;
  delete [] m0p;
}


void HelmholtzTable::gaussq(int N) {
  // First compute the nodes x_p and weights w_p for Gaussian quadrature to
  // approximate \int_{-1}^1 f(x) dx
  double muzero = 2.0;
  double *B = new double[N];
  double NM1 = N - 1;

  for (int i = 1; i <= NM1; ++i)
    B[i - 1] = i / sqrt(4 * i * i - 1);

  assert(imtql2(N, x_p_, B, w_p_) == 0);

  for (int i = 0; i < N; ++i)
    w_p_[i] = muzero * w_p_[i] * w_p_[i];

  // Scale to approximate \int_0^{\pi / 2} f(x) dx
  for (int i = 0; i < N; ++i) {
    x_p_[i] = (x_p_[i] + 1) / 4 * M_PI;
    w_p_[i] = w_p_[i] * sin(x_p_[i]) * M_PI / 4;
  }

  delete [] B;
}

HelmholtzTable::~HelmholtzTable() {
  delete [] sqf_;
  for (auto it = dmat_plus_->begin(); it != dmat_plus_->end(); ++it) {
    delete [] it->second;
  }
  for (auto it = dmat_minus_->begin(); it != dmat_minus_->end(); ++it) {
    delete [] it->second;
  }
  delete dmat_plus_;
  delete dmat_minus_;
  delete [] m2m_;
  delete [] l2l_;
  delete [] x_e_;
  delete [] w_e_;
  delete [] m_e_;
  delete [] sm_e_;
  delete [] n_e_;
  delete [] f_e_;
  delete [] smf_e_;
  delete [] x_p_;
  delete [] w_p_;
  delete [] m_p_;
  delete [] sm_p_;
  delete [] n_p_;
  delete [] f_p_;
  delete [] smf_p_;

  for (int i = 0; i <= maxlev; ++i) {
    delete [] ealphaj_e_[i];
    delete [] xs_e_[i];
    delete [] ys_e_[i];
    delete [] zs_e_[i];
    delete [] ealphaj_p_[i];
    delete [] xs_p_[i];
    delete [] ys_p_[i];
    delete [] zs_p_[i];
  }

  delete [] ealphaj_e_;
  delete [] xs_e_;
  delete [] ys_e_;
  delete [] zs_e_;
  delete [] ealphaj_p_;
  delete [] xs_p_;
  delete [] ys_p_;
  delete [] zs_p_;
}

int HelmholtzTable::imtql2(int N, double *D, double *E, double *Z) {
  // This function computes the eigenvalues and first component of the
  // eigenvector of a symmetric tridiagonal matrix by the implicit QL method.
  // On input, N is the order of the matrix, d has the diagonal elements, and e
  // has the subdigonal elements of the matrix in its first N - 1 position.
  // On output, d contains the eigenvalues in ascending order, and z contains
  // the first components of the orthonormal eigenvectors of the symmetric
  // tridiagonal matrix.
  // The return value is 0 for a normal return and j if the j-th eigenvalue has
  // not been determined after 30 iterations.

  // Initialize z to be the first row of the identity matrix
  Z[0] = 1;

  const int max_iter = 30;
  const double epsilon = 1e-14;
  int err = 0;

  if (N == 1)
    return err;

  for (int L = 1; L <= N; ++L) {
    int J = 0;

    while (105) {
      // Look for small sub-diagonal element
      int M;
      for (M = L; M <= N; ++M) {
        if (M == N)
          break;

        if (fabs(E[M - 1]) <= epsilon * (fabs(D[M - 1]) + fabs(D[M])))
          break;
      }

      double P = D[L - 1];

      if (M == L)
        break;

      if (J == max_iter) {
        err = J;
        return err;
      }

      J++;

      // Form shift
      double G = (D[L] - P) / (2 * E[L - 1]);
      double R = sqrt(G * G + 1);
      if (G >= 0) {
        G = D[M - 1] - P + E[L - 1] / (G + fabs(R));
      } else {
        G = D[M - 1] - P + E[L - 1] / (G - fabs(R));
      }
      double S = 1.0;
      double C = 1.0;
      P = 0.0;

      int MML = M - L;
      for (int II = 1; II <= MML; ++II) {
        int I = M - II;
        double F = S * E[I - 1];
        double B = C * E[I - 1];

        if (fabs(F) < fabs(G)) {
          S = F / G;
          R = sqrt(S * S + 1);
          E[I] = G * R;
          C = 1.0 / R;
          S = S * C;
        } else {
          C = G / F;
          R = sqrt(C * C + 1);
          E[I] = F * R;
          S = 1.0 / R;
          C = C * S;
        }

        G = D[I] - P;
        R = (D[I - 1] - G) * S + 2 * C * B;
        P = S * R;
        D[I] = G + P;
        G = C * R - B;

        // Form first component of vector
        F = Z[I];
        Z[I] = S * Z[I - 1] + C * F;
        Z[I - 1] = C * Z[I - 1] - S * F;
      }

      D[L - 1] -= P;
      E[L - 1] = G;
      E[M] = 0.0;
    }
  }

  for (int II = 2; II <= N; ++II) {
    int I = II - 1;
    int K = I;
    double P = D[I - 1];

    for (int j = II; j <= N; ++j) {
      if (D[j - 1] >= P)
        continue;
      K = j;
      P = D[j - 1];
    }

    if (K == I)
      continue;

    D[K - 1] = D[I - 1];
    D[I - 1] = P;
    P = Z[I - 1];
    Z[I - 1] = Z[K - 1];
    Z[K - 1] = P;
  }

  return err;
}

void HelmholtzTable::generate_sqf() {
  sqf_ = new double[(p_ + 1) * (p_ + 2) / 2];
  double *temp = new double[2 * p_ + 1];

  temp[0] = 1.0;
  for (int i = 1; i <= p_ * 2; ++i) {
    temp[i] = temp[i - 1] * i;
  }

  int idx = 0;
  for (int n = 0; n <= p_; ++n) {
    int n1 = 2 * n + 1;
    for (int m = 0; m <= n; ++m) {
      sqf_[idx++] = n1 * temp[n - m] / temp[n + m];
    }
  }

  delete [] temp;
}

void HelmholtzTable::generate_scaled_wigner_dmat() {
  dmat_plus_ = new builtin_map_t;
  dmat_minus_ = new builtin_map_t;
  double cbeta[3] = {sqrt(3) / 3, -sqrt(3) / 3, 0};
  //int nd = (p_ + 1) * (4 * p_ * p_ + 11 * p_ + 6) / 6;
  int nd = (p_ + 1) * (2 * p_ + 1) * (2 * p_ + 3) / 3; 
  for (int i = 0; i < 3; ++i) {
    double beta = acos(cbeta[i]);
    double *dp_data = new double[nd];
    double *dm_data = new double[nd];
    generate_scaled_dmat_of_beta(beta, dp_data, dm_data);
    (*dmat_plus_)[cbeta[i]] = dp_data;
    (*dmat_minus_)[cbeta[i]] = dm_data;
  }
}

void HelmholtzTable::generate_scaled_dmat_of_beta(double beta, double *dp,
                                                  double *dm) {
  double cbeta = cos(beta);
  double sbeta = sin(beta);
  double s2beta2 = (1 - cbeta) / 2; // sin^2(beta / 2)
  double c2beta2 = (1 + cbeta) / 2; // cos^2(beta / 2)

  // Set d_0^{0, 0} to 1
  dp[0] = 1;
  dm[0] = 1;

  // Set d_1^{0, m}
  dp[hidx(1, 0, -1)] = -sbeta / sqrt(2); // d_1^{0, -1}
  dp[hidx(1, 0, 0)] = cbeta; // d_1^{0, 0}
  dp[hidx(1, 0, 1)] = sbeta / sqrt(2); // d_1^{0, 1}
  dm[hidx(1, 0, -1)] = -dp[hidx(1, 0, -1)];
  dm[hidx(1, 0, 0)] = dp[hidx(1, 0, 0)];
  dm[hidx(1, 0, 1)] = -dp[hidx(1, 0, 1)];

  // Set d_1^{1, m}
  dp[hidx(1, 1, -1)] = s2beta2; // d_1^{1, -1}
  dp[hidx(1, 1, 0)] = -sbeta / sqrt(2); // d_1^{1, 0}
  dp[hidx(1, 1, 1)] = c2beta2; // d_1^{1, 1}
  dm[hidx(1, 1, -1)] = dp[hidx(1, 1, -1)];
  dm[hidx(1, 1, 0)] = -dp[hidx(1, 1, 0)];
  dm[hidx(1, 1, 1)] = dp[hidx(1, 1, 1)];


  // Compute d_n^{0, m} for 2 <= n <= P
  for (int n = 2; n <= p_; ++n) {
    double *dpc = NULL, *dpp = NULL, *dmc = NULL;
    int m;
    // Get address of d_n^{0, 0} for angle beta, saved in dpc
    // Get address of d_{n - 1}^{0, 0} for angle beta, saved in dpp
    // Get address of d_n^{0, 0} for angle -beta, saved in dmc
    dpc = &dp[hidx(n, 0, 0)];
    dpp = &dp[hidx(n - 1, 0, 0)];
    dmc = &dm[hidx(n, 0, 0)];


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

      dpc = &dp[hidx(n, mp, 0)];
      dpp = &dp[hidx(n - 1, mp - 1, 0)];
      dmc = &dm[hidx(n, mp, 0)];

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

  // Fill d_n^{mp, m} for negative mp values 
  for (int n = 1; n <= p_; ++n) {
    for (int mp = -n; mp <= -1; ++mp) {
      for (int m = -n; m <= -1; ++m) {
        dp[hidx(n, mp, m)] = dp[hidx(n, -m, -mp)]; 
        dm[hidx(n, mp, m)] = dm[hidx(n, -m, -mp)]; 
      } 

      for (int m = 0; m <= n; ++m) {
        dp[hidx(n, mp, m)] = pow(-1, m - mp) * dp[hidx(n, m, mp)]; 
        dm[hidx(n, mp, m)] = pow(-1, m - mp) * dm[hidx(n, m, mp)]; 
      }
    }
  }

  // Scale
  double factorial[2 * p_ + 1];
  factorial[0] = 1.0;
  for (int i = 1; i <= 2 * p_; ++i) {
    factorial[i] = factorial[i - 1] * i;
  }

  int offset = 0; 
  for (int n = 0; n <= p_; ++n) {
    for (int mp = -n; mp <= n; ++mp) {
      for (int m = -n; m <= n; ++m) {
        double scale = sqrt(factorial[n - abs(mp)] * factorial[n + abs(m)] / 
                            factorial[n + abs(mp)] / factorial[n - abs(m)]); 
        dp[offset] *= scale; 
        dm[offset] *= scale; 
        offset++;
      }
    }
  }
}

void HelmholtzTable::generate_m2m() {
  // Need to compute coefficients for boxes from levels 3 to maxlev
  m2m_ = new double[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * (maxlev - 2)];
  double *factorial = new double[2 * p_ + 1];
  double *bessel = new double[2 * p_ + 1];

  factorial[0] = 1.0;
  for (int i = 1; i <= 2 * p_; ++i) {
    factorial[i] = factorial[i - 1] * i;
  }

  for (int lev = 3; lev <= maxlev; ++lev) {
    // Compute shift distance rho
    double rho = sqrt(3) / 2 * size_ / pow(2, lev);

    // Compute scaling factor
    double sigma = scale_ / pow(2, lev);

    // Compute scaled spehrical bessel function of the first kind
    bessel_jn_scaled(2 * p_, omega_ * rho, omega_ * rho, bessel);

    double *coeff = &m2m_[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * (lev - 3)];

    for (int n = 0; n <= p_; ++n) {
      for (int m = 0; m <= n; ++m) {
        for (int np = m; np <= p_; ++np) {
          double temp = 0;
          int bound = (n <= np ? n : np);
          for (int k = m; k <= bound; ++k) {
            temp += pow(-0.5, n + k) * pow(sigma, np - n) *
              (2 * n + 1) * factorial[n - m] * factorial[np + m] *
              factorial[2 * k] / factorial[k + m] / factorial[k] /
              factorial[k - m] / factorial[n - k] / factorial[np - k] *
              bessel[np + n - k] * pow(omega_ * rho, np + n - 2 * k);
          }
          coeff[sidx(n, m, np, p_)] = temp;
        }
      }
    }
  }

  delete [] factorial;
  delete [] bessel;
}

void HelmholtzTable::generate_l2l() {
  // Need to compute coefficient for boxes from levels 2 to (maxlev - 1)
  l2l_ = new double[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * (maxlev - 2)];
  double *factorial = new double[2 * p_ + 1];
  double *bessel = new double[2 * p_ + 1];

  for (int lev = 2; lev <= maxlev - 1; ++lev) {
    // Compute shift distance rho
    double rho = sqrt(3) / 4 * size_ / pow(2, lev);

    // Compute scaling factor
    double sigma = scale_ / pow(2, lev);

    // Compute scaled spherical bessel function of the first kind
    bessel_jn_scaled(2 * p_, omega_ * rho, omega_ * rho, bessel);

    double *coeff = &l2l_[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * (lev - 2)];

    for (int n = 0; n <= p_; ++n) {
      for (int m = 0; m <= n; ++m) {
        for (int np = m; np <= p_; ++np) {
          double temp = 0;
          int bound = (n <= np ? n : np);
          for (int k = m; k <= bound; ++k) {
            temp += pow(-0.5, n + k) * pow(sigma, n - np) *
              (2 * n + 1) * factorial[n - m] / factorial[k + m] *
              factorial[np + m] * factorial[2 * k] / factorial[k] /
              factorial[k - m] / factorial[n - k] / factorial[np - k] *
              bessel[np + n - k] * pow(omega_ * rho, n + np - 2 * k);
          }
          coeff[sidx(n, m, np, p_)] = temp;
        }
      }
    }
  }


  delete [] factorial;
  delete [] bessel;
}

void update_helmholtz_table(int n_digits, double size, double omega) {
  if (builtin_helmholtz_table_ == nullptr)
    builtin_helmholtz_table_ = std::unique_ptr<HelmholtzTable>{
      new HelmholtzTable{n_digits, size, omega}};
}



} // namespace dashmm

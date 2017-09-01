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


/// \file
/// Implementation of Laplace kernel


#include "builtins/laplace.h"

namespace dashmm {

void lap_rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR) {
  int p = builtin_laplace_table_->p();
  
  // Compute exp(i * alpha)
  dcomplex_t ealpha{cos(alpha), sin(alpha)};
  
  // Compute powers of exp(i * alpha)
  std::vector<dcomplex_t>powers_ealpha(p + 1); 
  powers_ealpha[0] = dcomplex_t{1.0, 0.0};
  for (int j = 1; j <= p; ++j) 
    powers_ealpha[j] = powers_ealpha[j - 1] * ealpha;
  
  int offset = 0;
  for (int n = 0; n <= p; n++) {
    for (int m = 0; m <= n; m++) {
      MR[offset] = M[offset] * powers_ealpha[m];
      offset++;
    }
  }
}

void lap_rotate_sph_y(const dcomplex_t *M, const double *d, dcomplex_t *MR) {
  int p = builtin_laplace_table_->p();
  
  int offset = 0;
  for (int n = 0; n <= p; ++n) {
    int power_mp = 1;
    for (int mp = 0; mp <= n; ++mp) {
      // Retrieve address of wigner d-matrix entry d_n^{mp, 0}
      const double *coeff = &d[didx(n, mp, 0)];
      // Get address of original harmonic expansion M_n^0
      const dcomplex_t *Mn = &M[midx(n, 0)];
      // Compute rotated spherical harmonic M_n^mp
      MR[offset] = Mn[0] * coeff[0];
      double power_m = -1;
      for (int m = 1; m <= n; ++m) {
        MR[offset] += (Mn[m] * power_m * coeff[m] + conj(Mn[m]) * coeff[-m]);
        power_m = -power_m;
      }
      MR[offset++] *= power_mp;
      power_mp = -power_mp;
    }
  }  
}

void lap_s_to_m(Point dist, double q, double scale, dcomplex_t *M) {
  int p = builtin_laplace_table_->p();
  const double *sqf = builtin_laplace_table_->sqf();

  std::vector<double> legendre((p + 1) * (p + 2) / 2); 
  std::vector<double> powers_r(p + 1); 
  std::vector<dcomplex_t> powers_ephi(p + 1); 

  powers_r[0] = 1.0; 
  powers_ephi[0] = dcomplex_t{1.0, 0.0}; 

  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();
  double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

  // Compute exp(-i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, -dist.y() / proj});

  // Compute powers of r
  r *= scale;
  for (int j = 1; j <= p; ++j) 
    powers_r[j] = powers_r[j - 1] * r;
  
  // Compute powers of exp(-i * phi)
  for (int j = 1; j <= p; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;
  
  // Compute multipole expansion M_n^m
  legendre_Plm(p, ctheta, legendre.data());
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      M[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
        legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
    }
  }
}

void lap_s_to_l(Point dist, double q, double scale, dcomplex_t *L) {
  int p = builtin_laplace_table_->p();
  const double *sqf = builtin_laplace_table_->sqf();

  std::vector<double> legendre((p + 1) * (p + 2) / 2); 
  std::vector<double> powers_r(p + 1); 
  std::vector<dcomplex_t> powers_ephi(p + 1); 
  powers_ephi[0] = dcomplex_t{1.0, 0.0};

  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();

  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);
  
  // Compute exp(-i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, -dist.y() / proj});

  // Compute powers of 1 / r
  powers_r[0] = 1.0 / r;
  r *= scale;
  for (int j = 1; j <= p; ++j) 
    powers_r[j] = powers_r[j - 1] / r;
  
  // Compute powers of exp(-i * phi)
  for (int j = 1; j <= p; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;
  
  // compute local expansion L_n^m
  legendre_Plm(p, ctheta, legendre.data());
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      L[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
        legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
    }
  }
}

void lap_m_to_m(int from_child, const dcomplex_t *M, dcomplex_t *W) {
  int p = builtin_laplace_table_->p();
  const double *sqbinom = builtin_laplace_table_->sqbinom();

  // Get precomputed Wigner d-matrix for rotation about the y-axis
  const double *d1 = (from_child < 4 ?
                      builtin_laplace_table_->dmat_plus(1.0 / sqrt(3.0)) :
                      builtin_laplace_table_->dmat_plus(-1.0 / sqrt(3.0)));
  const double *d2 = (from_child < 4 ?
                      builtin_laplace_table_->dmat_minus(1.0 / sqrt(3.0)) :
                      builtin_laplace_table_->dmat_minus(-1.0 / sqrt(3.0)));
  
  // Shift distance along the z-axis, combined with Y_n^0(pi, 0)
  const double rho = -sqrt(3) / 2;

  // Compute powers of rho
  std::vector<double> powers_rho(p + 1); 
  powers_rho[0] = 1.0;
  for (int i = 1; i <= p; ++i) {
    powers_rho[i] = powers_rho[i - 1] * rho;
  }

  std::vector<dcomplex_t> T((p + 1) * (p + 2) / 2); 

  // Table of rotation angle about the z-axis, as an integer multiple of pi/4
  const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
  // Get rotation angle about the z-axis
  double alpha = tab_alpha[from_child] * M_PI_4;

  // Rotate the multipole expansion of the child box about z-axis
  lap_rotate_sph_z(M, alpha, W); 

  // Rotate the previous result further about the y-axis
  lap_rotate_sph_y(W, d1, T.data());

  // Offset to operate multipole expansion
  int offset = 0;

  // Shift along the z-axis by a distance of rho, write result in W
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      W[offset] = T[offset];
      for (int k = 1; k <= n - m; ++k) {
        W[offset] += T[midx(n - k, m)] * powers_rho[k] *
          sqbinom[midx(n - m, k)] * sqbinom[midx(n + m, k)];
      }
      offset++;
    }
  }
  
  // Reverse rotate the shifted harmonic expansion about the y-axis
  lap_rotate_sph_y(W, d2, T.data());

  // Reverse rotate the previous result further about the z-axis
  lap_rotate_sph_z(T.data(), -alpha, W);

  double temp = 1;
  offset = 0;
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      W[offset++] *= temp;
    }
    temp /= 2;
  }
}

void lap_l_to_l(int to_child, const dcomplex_t *L, dcomplex_t *W) {
  int p = builtin_laplace_table_->p();
  const double *sqbinom = builtin_laplace_table_->sqbinom();

  // Get precomputed Wigner d-matrix for rotation about the y-axis
  const double *d1 = (to_child < 4 ?
                      builtin_laplace_table_->dmat_plus(1.0 / sqrt(3.0)) :
                      builtin_laplace_table_->dmat_plus(-1.0 / sqrt(3.0)));
  const double *d2 = (to_child < 4 ?
                      builtin_laplace_table_->dmat_minus(1.0 / sqrt(3.0)) :
                      builtin_laplace_table_->dmat_minus(-1.0 / sqrt(3.0)));

  // Shift distance along the z-axis, combined with Y_n^0(pi, 0)
  const double rho = -sqrt(3) / 4;

  // Compute powers of rho
  std::vector<double> powers_rho(p + 1); 
  powers_rho[0] = 1.0;
  for (int i = 1; i <= p; ++i) {
    powers_rho[i] = powers_rho[i - 1] * rho;
  }

  std::vector<dcomplex_t> T((p + 1) * (p + 2) / 2);

  // Table of rotation angle about the z-axis as an integer multiple of pi / 4
  const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
  // Get rotation angle about the z-axis
  double alpha = tab_alpha[to_child] * M_PI_4;
  
  // Rotate the local expansion of the parent box about z-axis
  lap_rotate_sph_z(L, alpha, W);

  // Rotate the previous result further about the y-axis
  lap_rotate_sph_y(W, d1, T.data());

  // Offset to operate local expansion
  int offset = 0;
  
  // Shift along the z-axis by a distance of rho, write result in W1
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      W[offset] = T[offset];
      for (int k = 1; k <= p - n; k++) {
        W[offset] += T[midx(n + k, m)] * powers_rho[k] *
          sqbinom[midx(n + k - m, k)] * sqbinom[midx(n + k + m, k)];
      }
      offset++;
    }
  }
  
  // Reverse rotate the shifted harmonic expansion about the y-axis
  lap_rotate_sph_y(W, d2, T.data());
  
  // Reverse rotate the previous result further about the z-axis
  lap_rotate_sph_z(T.data(), -alpha, W);
  
  double temp = 1;
  offset = 0;
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      W[offset++] *= temp;
    }
    temp /= 2;
  }
}

void lap_m_to_i(const dcomplex_t *M, ViewSet &views, int id) {
  // Addresses of the views
  dcomplex_t *E_px = reinterpret_cast<dcomplex_t *>(views.view_data(id));
  dcomplex_t *E_mx = reinterpret_cast<dcomplex_t *>(views.view_data(id + 1));
  dcomplex_t *E_py = reinterpret_cast<dcomplex_t *>(views.view_data(id + 2));
  dcomplex_t *E_my = reinterpret_cast<dcomplex_t *>(views.view_data(id + 3));
  dcomplex_t *E_pz = reinterpret_cast<dcomplex_t *>(views.view_data(id + 4));
  dcomplex_t *E_mz = reinterpret_cast<dcomplex_t *>(views.view_data(id + 5));

  // Addresses of exponential expansions in the positive axis direction
  dcomplex_t *EP[3] = {E_px, E_py, E_pz};
  // Addresses of exponential expansions in the negative axis direction
  dcomplex_t *EM[3] = {E_mx, E_my, E_mz};

  int p = builtin_laplace_table_->p();
  int s = builtin_laplace_table_->s();
  int nsh = (p + 1) * (p + 2) / 2;
  const double *weight_ = builtin_laplace_table_->weight();
  const int *m_ = builtin_laplace_table_->m();
  const int *f_ = builtin_laplace_table_->f();
  const int *smf_ = builtin_laplace_table_->smf();
  const double *d1 = builtin_laplace_table_->dmat_plus(0.0);
  const double *d2 = builtin_laplace_table_->dmat_minus(0.0);
  const double *lambdaknm = builtin_laplace_table_->lambdaknm();
  const dcomplex_t *ealphaj = builtin_laplace_table_->ealphaj();

  // Allocate scratch space to handle x- and y-direction
  // exponential expansions
  std::vector<dcomplex_t> W1((p + 1) * (p + 2) / 2);
  std::vector<dcomplex_t> W2((p + 1) * (p + 2) / 2);
  
  // Setup y-direction. Rotate the multipole expansion M about z-axis by -pi /
  // 2, making (x, y, z) frame (-y, x, z). Next, rotate it again about the new
  // y axis by -pi / 2. The (-y, x, z) in the first rotated frame becomes (z,
  // x, y) in the final frame.
  lap_rotate_sph_z(M, -M_PI / 2, W1.data());
  lap_rotate_sph_y(W1.data(), d2, W2.data());

  // Setup x-direction. Rotate the multipole expansion M about y axis by pi /
  // 2. This makes (x, y, z) frame into (-z, y, x).
  lap_rotate_sph_y(M, d1, W1.data());

  // Addresses of the spherical harmonic expansions
  const dcomplex_t *SH[3] = {W1.data(), W2.data(), M};

  for (int dir = 0; dir <= 2; ++dir) {
    int offset = 0;
    for (int k = 0; k < s ; ++k) {
      double weight = weight_[k] / m_[k];

      // Compute sum_{n = m}^p M_n^m * lambda_k^n / sqrt((n+m)! * (n - m)!)
      // z1 handles M_n^m where n is even
      std::vector<dcomplex_t> z1(f_[k] + 1);
      // z2 handles M_n^m where n is odd
      std::vector<dcomplex_t> z2(f_[k] + 1);
      
      // Process M_n^0 terms
      z1[0] = 0;
      z2[0] = 0;
      for (int n = 0; n <= p; n += 2) {
        z1[0] += SH[dir][midx(n, 0)] * lambdaknm[k * nsh + midx(n, 0)];
      }
      for (int n = 1; n <= p; n += 2) {
        z2[0] += SH[dir][midx(n, 0)] * lambdaknm[k * nsh + midx(n, 0)];
      }

      // Process M_n^m terms for nonzero m
      for (int m = 1; m <= f_[k]; m += 2) {
        z1[m] = 0;
        z2[m] = 0;
        for (int n = m; n <= p; n += 2) {
          z2[m] += SH[dir][midx(n, m)] * lambdaknm[k * nsh + midx(n, m)];
        }
        for (int n = m + 1; n <= p; n += 2) {
          z1[m] += SH[dir][midx(n, m)] * lambdaknm[k * nsh + midx(n, m)];
        }
      }
      
      for (int m = 2; m <= f_[k]; m += 2) {
        z1[m] = 0;
        z2[m] = 0;
        for (int n = m; n <= p; n += 2) {
          z1[m] += SH[dir][midx(n, m)] * lambdaknm[k * nsh + midx(n, m)];
        }
        for (int n = m + 1; n <= p; n += 2) {
          z2[m] += SH[dir][midx(n, m)] * lambdaknm[k * nsh + midx(n, m)];
        }
      }
      
      // Compute W(k, j)
      for (int j = 1; j <= m_[k] / 2; ++j) {
        dcomplex_t up = z1[0] + z2[0]; // accumulate +dir
        dcomplex_t dn = z1[0] - z2[0]; // accumulate -dir
        dcomplex_t power_I {0.0, 1.0};
        for (int m = 1; m <= f_[k]; ++m) {
          int idx = smf_[k] + (j - 1) * f_[k] + (m - 1);
          up += 2 * real(ealphaj[idx] * (z1[m] + z2[m])) * power_I;
          dn += 2 * real(ealphaj[idx] * (z1[m] - z2[m])) * power_I;
          power_I *= dcomplex_t{0.0, 1.0};
        }
        EP[dir][offset] = weight * up;
        EM[dir][offset] = weight * dn;
        offset++;
      }
    }
  }
}

void lap_i_to_i(Index s_index, Index t_index, const ViewSet &s_views, 
                int sid, int tid, ViewSet &t_views) {
  const dcomplex_t *S[6]{
    reinterpret_cast<dcomplex_t *>(s_views.view_data(sid)),
      reinterpret_cast<dcomplex_t *>(s_views.view_data(sid + 1)),
      reinterpret_cast<dcomplex_t *>(s_views.view_data(sid + 2)),
      reinterpret_cast<dcomplex_t *>(s_views.view_data(sid + 3)),
      reinterpret_cast<dcomplex_t *>(s_views.view_data(sid + 4)),
      reinterpret_cast<dcomplex_t *>(s_views.view_data(sid + 5))
      };
  
  // Compute index offsets between the current source node and the 1st child
  // of the parent node
  int dx = s_index.x() - t_index.x() * 2;
  int dy = s_index.y() - t_index.y() * 2;
  int dz = s_index.z() - t_index.z() * 2;
  
  // Exponential expansions on the source side
  int nexp = builtin_laplace_table_->nexp();

  // Each S is going to generate between 1 and 3 views of the exponential
  // expansions on the target side.
  size_t view_size = nexp * sizeof(dcomplex_t);
  dcomplex_t *T1 = new dcomplex_t[nexp]();
  dcomplex_t *T2 = new dcomplex_t[nexp]();
  dcomplex_t *T3 = new dcomplex_t[nexp]();
  char *C1 = reinterpret_cast<char *>(T1);
  char *C2 = reinterpret_cast<char *>(T2);
  char *C3 = reinterpret_cast<char *>(T3);
  dcomplex_t *T[3] = {T1, T2, T3};
  char *C[3] = {C1, C2, C3};
  bool used[3] = {false, false, false};

  for (int i = 0; i < 3; ++i) {
    int tag = merge_and_shift_table[dx + 2][dy + 2][dz + 2][i];

    if (tag == -1)
      break;

    if (tag <= 1) {
      lap_e_to_e(T[i], S[5], dx, dy, 0);
    } else if (tag <= 5) {
      lap_e_to_e(T[i], S[3], dz, dx, 0);
    } else if (tag <= 13) {
      lap_e_to_e(T[i], S[1], -dz, dy, 0);
    } else if (tag <= 15) {
      lap_e_to_e(T[i], S[4], -dx, -dy, 0);
    } else if (tag <= 19) {
      lap_e_to_e(T[i], S[2], -dz, -dx, 0);
    } else {
      lap_e_to_e(T[i], S[0], dz, -dy, 0);
    }
    
    t_views.add_view(tid + tag, view_size, C[i]);
    used[i] = true;
  }
  
  if (used[1] == false)
    delete [] T2;
  
  if (used[2] == false)
    delete [] T3;
}

void lap_i_to_l(const ViewSet &views, int id, Index t_index, 
                double scale, dcomplex_t *L) {
  const dcomplex_t *E[28]{
    reinterpret_cast<dcomplex_t *>(views.view_data(id)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 1)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 2)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 3)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 4)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 5)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 6)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 7)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 8)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 9)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 10)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 11)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 12)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 13)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 14)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 15)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 16)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 17)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 18)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 19)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 20)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 21)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 22)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 23)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 24)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 25)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 26)),
      reinterpret_cast<dcomplex_t *>(views.view_data(id + 27))
      };

  int to_child = 4 * (t_index.z() % 2) + 2 * (t_index.y() % 2) +
    (t_index.x() % 2);
  
  int nexp = builtin_laplace_table_->nexp();
  dcomplex_t *S = new dcomplex_t[nexp * 6]();
  dcomplex_t *S_mz = S;
  dcomplex_t *S_pz = S + nexp;
  dcomplex_t *S_my = S + 2 * nexp;
  dcomplex_t *S_py = S + 3 * nexp;
  dcomplex_t *S_mx = S + 4 * nexp;
  dcomplex_t *S_px = S + 5 * nexp;
  
  switch (to_child) {
  case 0:
    lap_e_to_e(S_mz, E[uall], 0, 0, 3);
    lap_e_to_e(S_mz, E[u1234], 0, 0, 2);
    lap_e_to_e(S_pz, E[dall], 0, 0, 2);
    
    lap_e_to_e(S_my, E[nall], 0, 0, 3);
    lap_e_to_e(S_my, E[n1256], 0, 0, 2);
    lap_e_to_e(S_my, E[n12], 0, 0, 2);
    lap_e_to_e(S_py, E[sall], 0, 0, 2);
    
    lap_e_to_e(S_mx, E[eall], 0, 0, 3);
    lap_e_to_e(S_mx, E[e1357], 0, 0, 2);
    lap_e_to_e(S_mx, E[e13], 0, 0, 2);
    lap_e_to_e(S_mx, E[e1], 0, 0, 2);
    lap_e_to_e(S_px, E[wall], 0, 0, 2);
    break;
  case 1:
    lap_e_to_e(S_mz, E[uall], -1, 0, 3);
    lap_e_to_e(S_mz, E[u1234], -1, 0, 2);
    lap_e_to_e(S_pz, E[dall], 1, 0, 2);
    
    lap_e_to_e(S_my, E[nall], 0, -1, 3);
    lap_e_to_e(S_my, E[n1256], 0, -1, 2);
    lap_e_to_e(S_my, E[n12], 0, -1, 2);
    lap_e_to_e(S_py, E[sall], 0, 1, 2);
    
    lap_e_to_e(S_mx, E[eall], 0, 0, 2);
    lap_e_to_e(S_px, E[wall], 0, 0, 3);
    lap_e_to_e(S_px, E[w2468], 0, 0, 2);
    lap_e_to_e(S_px, E[w24], 0, 0, 2);
    lap_e_to_e(S_px, E[w2], 0, 0, 2);
    break;
  case 2:
    lap_e_to_e(S_mz, E[uall], 0, -1, 3);
    lap_e_to_e(S_mz, E[u1234], 0, -1, 2);
    lap_e_to_e(S_pz, E[dall], 0, 1, 2);
    
    lap_e_to_e(S_my, E[nall], 0, 0, 2);
    lap_e_to_e(S_py, E[sall], 0, 0, 3);
    lap_e_to_e(S_py, E[s3478], 0, 0, 2);
    lap_e_to_e(S_py, E[s34], 0, 0, 2);
    
    lap_e_to_e(S_mx, E[eall], 0, -1, 3);
    lap_e_to_e(S_mx, E[e1357], 0, -1, 2);
    lap_e_to_e(S_mx, E[e13], 0, -1, 2);
    lap_e_to_e(S_mx, E[e3], 0, -1, 2);
    lap_e_to_e(S_px, E[wall], 0, 1, 2);
    break;
  case 3:
    lap_e_to_e(S_mz, E[uall], -1, -1, 3);
    lap_e_to_e(S_mz, E[u1234], -1, -1, 2);
    lap_e_to_e(S_pz, E[dall], 1, 1, 2);
    
    lap_e_to_e(S_my, E[nall], 0, -1, 2);
    lap_e_to_e(S_py, E[sall], 0, 1, 3);
    lap_e_to_e(S_py, E[s3478], 0, 1, 2);
    lap_e_to_e(S_py, E[s34], 0, 1, 2);
    
    lap_e_to_e(S_mx, E[eall], 0, -1, 2);
    lap_e_to_e(S_px, E[wall], 0, 1, 3);
    lap_e_to_e(S_px, E[w2468], 0, 1, 2);
    lap_e_to_e(S_px, E[w24], 0, 1, 2);
    lap_e_to_e(S_px, E[w4], 0, 1, 2);
    break;
  case 4:
    lap_e_to_e(S_mz, E[uall], 0, 0, 2);
    lap_e_to_e(S_pz, E[dall], 0, 0, 3);
    lap_e_to_e(S_pz, E[d5678], 0, 0, 2);
    
    lap_e_to_e(S_my, E[nall], -1, 0, 3);
    lap_e_to_e(S_my, E[n1256], -1, 0, 2);
    lap_e_to_e(S_my, E[n56], -1, 0, 2);
    lap_e_to_e(S_py, E[sall], 1, 0, 2);
    
    lap_e_to_e(S_mx, E[eall], 1, 0, 3);
    lap_e_to_e(S_mx, E[e1357], 1, 0, 2);
    lap_e_to_e(S_mx, E[e57], 1, 0, 2);
    lap_e_to_e(S_mx, E[e5], 1, 0, 2);
    lap_e_to_e(S_px, E[wall], -1, 0, 2);
    break;
  case 5:
    lap_e_to_e(S_mz, E[uall], -1, 0, 2);
    lap_e_to_e(S_pz, E[dall], 1, 0, 3);
    lap_e_to_e(S_pz, E[d5678], 1, 0, 2);
    
    lap_e_to_e(S_my, E[nall], -1, -1, 3);
    lap_e_to_e(S_my, E[n1256], -1, -1, 2);
    lap_e_to_e(S_my, E[n56], -1, -1, 2);
    lap_e_to_e(S_py, E[sall], 1, 1, 2);
    
    lap_e_to_e(S_mx, E[eall], 1, 0, 2);
    lap_e_to_e(S_px, E[wall], -1, 0, 3);
    lap_e_to_e(S_px, E[w2468], -1, 0, 2);
    lap_e_to_e(S_px, E[w68], -1, 0, 2);
    lap_e_to_e(S_px, E[w6], -1, 0, 2);
    break;
  case 6:
    lap_e_to_e(S_mz, E[uall], 0, -1, 2);
    lap_e_to_e(S_pz, E[dall], 0, 1, 3);
    lap_e_to_e(S_pz, E[d5678], 0, 1, 2);
    
    lap_e_to_e(S_my, E[nall], -1, 0, 2);
    lap_e_to_e(S_py, E[sall], 1, 0, 3);
    lap_e_to_e(S_py, E[s3478], 1, 0, 2);
    lap_e_to_e(S_py, E[s78], 1, 0, 2);
    
    lap_e_to_e(S_mx, E[eall], 1, -1, 3);
    lap_e_to_e(S_mx, E[e1357], 1, -1, 2);
    lap_e_to_e(S_mx, E[e57], 1, -1, 2);
    lap_e_to_e(S_mx, E[e7], 1, -1, 2);
    lap_e_to_e(S_px, E[wall], -1, 1, 2);
    break;
  case 7:
    lap_e_to_e(S_mz, E[uall], -1, -1, 2);
    lap_e_to_e(S_pz, E[dall], 1, 1, 3);
    lap_e_to_e(S_pz, E[d5678], 1, 1, 2);

    lap_e_to_e(S_my, E[nall], -1, -1, 2);
    lap_e_to_e(S_py, E[sall], 1, 1, 3);
    lap_e_to_e(S_py, E[s3478], 1, 1, 2);
    lap_e_to_e(S_py, E[s78], 1, 1, 2);
    
    lap_e_to_e(S_mx, E[eall], 1, -1, 2);
    lap_e_to_e(S_px, E[wall], -1, 1, 3);
    lap_e_to_e(S_px, E[w2468], -1, 1, 2);
    lap_e_to_e(S_px, E[w68], -1, 1, 2);
    lap_e_to_e(S_px, E[w8], -1, 1, 2);
    break;
  }

  for (int i = 0; i < 6 * nexp; ++i) {
    S[i] *= scale;
  }
  
  lap_e_to_l(S_mz, 'z', false, L);
  lap_e_to_l(S_pz, 'z', true, L);
  lap_e_to_l(S_my, 'y', false, L);
  lap_e_to_l(S_py, 'y', true, L);
  lap_e_to_l(S_mx, 'x', false, L);
  lap_e_to_l(S_px, 'x', true, L);
  
  delete [] S;
}

void lap_e_to_e(dcomplex_t *M, const dcomplex_t *W, int x, int y, int z) {
  const dcomplex_t *xs = builtin_laplace_table_->xs();
  const dcomplex_t *ys = builtin_laplace_table_->ys();
  const double *zs = builtin_laplace_table_->zs();
  int s = builtin_laplace_table_->s();
  const int *m = builtin_laplace_table_->m();
  const int *sm = builtin_laplace_table_->sm();
  
  int offset = 0;
  for (int k = 0; k < s; ++k) {
    // Shifting factor in z direction
    double factor_z = zs[4 * k + z];
    for (int j = 0; j < m[k] / 2; ++j) {
      int sidx = (sm[k] + j) * 7 + 3;
      // Shifting factor in x direction
      dcomplex_t factor_x = xs[sidx + x];
      // Shifting factor in y direction
      dcomplex_t factor_y = ys[sidx + y];
      M[offset] += W[offset] * factor_z * factor_y * factor_x;
      offset++;
    }
  }
}

void lap_e_to_l(const dcomplex_t *E, char dir, bool sgn, dcomplex_t *L) {
  const double *sqf = builtin_laplace_table_->sqf();
  const dcomplex_t *ealphaj = builtin_laplace_table_->ealphaj();
  const double *lambda = builtin_laplace_table_->lambda();
  const int *m_ = builtin_laplace_table_->m();
  const int *sm_ = builtin_laplace_table_->sm();
  const int *f_ = builtin_laplace_table_->f();
  const int *smf_ = builtin_laplace_table_->smf();
  int p = builtin_laplace_table_->p();
  int s = builtin_laplace_table_->s();

  dcomplex_t *contrib = nullptr;
  dcomplex_t *W1 = new dcomplex_t[(p + 1) * (p + 2) / 2];
  dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];
  
  for (int k = 0; k < s; ++k) {
    int Mk2 = m_[k] / 2;
    
    // Compute sum_{j = 1}^{M(k)} W(k, j) exp(-i * m * alpha_j)
    dcomplex_t *z = new dcomplex_t[f_[k] + 1];
    
    // m = 0
    z[0] = 0.0;
    for (int j = 1; j <= Mk2; ++j) {
      int idx = sm_[k] + j - 1;
      z[0] += (E[idx] + conj(E[idx]));
    }
    
    // m = 1, ..., F(k), where m is odd
    for (int m = 1; m <= f_[k]; m += 2) {
      z[m] = 0.0;
      for (int j = 1; j <= Mk2; ++j) {
        int idx1 = sm_[k] + j - 1;
        int idx2 = smf_[k] + (j - 1) * f_[k] + m - 1;
        z[m] += (E[idx1] - conj(E[idx1])) * conj(ealphaj[idx2]);
      }
    }
    
    // m = 2, ..., F(k), where m is even
    for (int m = 2; m <= f_[k]; m += 2) {
      z[m] = 0.0;
      for (int j = 1; j <= Mk2; ++j) {
        int idx1 = sm_[k] + j - 1;
        int idx2 = smf_[k] + (j - 1) * f_[k] + m - 1;
        z[m] += (E[idx1] + conj(E[idx1])) * conj(ealphaj[idx2]);
      }
    }

    // Compute lambda_k's contribution
    double power_lambdak = 1.0; // (-lambda_k)^n
    for (int n = 0; n <= p; ++n) {
      int mmax = fmin(n, f_[k]);
      for (int m = 0; m <= mmax; ++m) {
        W1[midx(n, m)] += power_lambdak * z[m];
      }
      power_lambdak *= -lambda[k];
    }
    delete [] z;
  }
  
  if (!sgn) {
    // If the exponential expansion is not along the positive direction of the
    // axis with respect to the source, flip the sign of the converted L_n^m
    // terms where n is odd.
    int offset = 1; // address of L_1^0
    for (int n = 1; n <= p; n += 2) {
      for (int m = 0; m <= n; ++m) {
        W1[offset++] *= -1;
      }
      offset += (n + 2); // skip the storage for the even value of n
    }
  }
  
  // Scale the local expansion by i^m / sqrt((n - m)!(n + m)!)
  int offset = 0;
  for (int n = 0; n <= p; ++n) {
    dcomplex_t power_I{1.0, 0.0};
    for (int m = 0; m <= n; ++m) {
      W1[offset++] *= power_I / sqf[n - m] / sqf[n + m];
      power_I *= dcomplex_t{0.0, 1.0};
    }
  }
  
  if (dir == 'z') {
    contrib = W1;
  } else if (dir == 'y') {
    const double *d = builtin_laplace_table_->dmat_plus(0.0);
    lap_rotate_sph_y(W1, d, W2);
    lap_rotate_sph_z(W2, M_PI / 2, W1);
    contrib = W1;
  } else if (dir == 'x') {
    const double *d = builtin_laplace_table_->dmat_minus(0.0);
    lap_rotate_sph_y(W1, d, W2);
    contrib = W2;
  }
  
  // Merge converted local expansion with the stored one
  offset = 0;
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      L[offset] += contrib[offset];
      offset++;
    }
  }
  
  delete [] W1;
  delete [] W2;
}

std::vector<double> lap_m_to_t(Point dist, double scale, 
                               const dcomplex_t *M, bool g) {
  std::vector<double> retval; 

  int p = builtin_laplace_table_->p();
  const double *sqf = builtin_laplace_table_->sqf();
  std::vector<double> legendre((p + 2) * (p + 3) / 2); 
  std::vector<double> powers_r(p + 2); 
  std::vector<dcomplex_t> powers_ephi(p + 2); 
  powers_ephi[0] = dcomplex_t{1.0, 0.0}; 

  // Compute potential first 
  dcomplex_t potential{0.0, 0.0}; 
  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();
  
  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

  // Compute exp(i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, dist.y() / proj});

  // Compute powers of 1 / r
  powers_r[0] = 1.0 / r;
  r *= scale;
  for (int j = 1; j <= p + 1; ++j) 
    powers_r[j] = powers_r[j - 1] / r;
  
  // Compute powers of exp(i * phi)
  for (int j = 1; j <= p + 1; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;

  // Evaluate the multipole expansion M_n^0
  legendre_Plm(p + 1, ctheta, legendre.data());
  for (int n = 0; n <= p; ++n) 
    potential += M[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)];


  // Evaluate the multipole expansions M_n^m, where m = 1, ..., p
  for (int n = 1; n <= p; ++n) {
    for (int m = 1; m <= n; ++m) {
      potential += 2.0 * real(M[midx(n, m)] * powers_ephi[m]) *
        powers_r[n] * legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
    }
  }

  retval.push_back(real(potential)); 

  if (g) {
    // Compute field values
    dcomplex_t zs1{0.0, 0.0}, zs2{0.0, 0.0}, zs3{0.0, 0.0}; 
    double fx = 0, fy = 0, fz = 0; 
    
    // M_0^0 term
    zs1 += powers_ephi[1] * real(M[midx(0, 0)]) * powers_r[1] * 
      legendre[midx(1, 1)];
    fz = real(M[midx(0, 0)]) * powers_r[1] * legendre[midx(1, 0)];

    // M_n^0 terms, n = 1, ..., p
    for (int n = 1; n <= p; ++n) {
      zs1 += powers_ephi[1] * real(M[midx(n, 0)]) * powers_r[n + 1] * 
        legendre[midx(n + 1, 1)]; 
      zs2 += M[midx(n, 1)] * powers_r[n + 1] * legendre[midx(n + 1, 0)] 
        * sqf[n + 1] / sqf[n - 1]; 
      fz += real(M[midx(n, 0)]) * powers_r[n + 1] * legendre[midx(n + 1, 0)]
        * (n + 1);
    }

    // M_n^m terms, n = 1, ..., p, m = 1, ..., n
    for (int n = 1; n <= p; ++n) {
      for (int m = 1; m <= n; ++m) {
        zs1 += M[midx(n, m)] * powers_r[n + 1] * powers_ephi[m + 1] * 
          legendre[midx(n + 1, m + 1)] * sqf[n - m] / sqf[n + m];
        if (m > 1) {
          zs2 += M[midx(n, m)] * powers_r[n + 1] * powers_ephi[m - 1] * 
            legendre[midx(n + 1, m - 1)] * sqf[n - m + 2] / sqf[n - m] *
            sqf[n - m + 2] / sqf[n + m]; 
        }
        zs3 += M[midx(n, m)] * powers_r[n + 1] * powers_ephi[m] * 
          legendre[midx(n + 1, m)] * sqf[n - m + 1] / sqf[n - m] * 
          sqf[n - m + 1] / sqf[n + m];         
      }
    }

    fx = real(zs2 - zs1); 
    fy = -imag(zs2 + zs1); 
    fz += 2.0 * real(zs3); 

    retval.push_back(fx / scale); 
    retval.push_back(fy / scale);
    retval.push_back(fz / scale); 
  }

  return retval; 
}

std::vector<double> lap_l_to_t(Point dist, double scale, 
                               const dcomplex_t *L, bool g) {
  std::vector<double> retval; 

  int p = builtin_laplace_table_->p();
  const double *sqf = builtin_laplace_table_->sqf();
  std::vector<double> legendre((p + 2) * (p + 3) / 2); 
  std::vector<double> powers_r(p + 2); 
  std::vector<dcomplex_t> powers_ephi(p + 2); 
  powers_r[0] = 1.0;
  powers_ephi[0] = dcomplex_t{1.0, 0.0};

  // Compute potential first
  dcomplex_t potential{0.0, 0.0};
  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();

  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

  // Compute exp(i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, dist.y() / proj});

  // Compute powers of r
  r *= scale;
  for (int j = 1; j <= p + 1; ++j) 
    powers_r[j] = powers_r[j - 1] * r;
  
  // Compute powers of exp(i * phi)
  for (int j = 1; j <= p + 1; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;
  
  // Evaluate the local expansion L_n^0
  legendre_Plm(p + 1, ctheta, legendre.data());
  for (int n = 0; n <= p; ++n) 
    potential += L[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)];
  
  // Evaluate the local expansions L_n^m, where m = 1, ..., p
  for (int n = 1; n <= p; ++n) {
    for (int m = 1; m <= n; ++m) {
      potential += 2.0 * real(L[midx(n, m)] * powers_ephi[m]) *
        powers_r[n] * legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
    }
  }

  retval.push_back(real(potential)); 

  if (g) {
    // Compute field values 
    dcomplex_t zs1{0.0, 0.0}, zs2{0.0, 0.0}, zs3{0.0, 0.0}; 
    double fx = 0, fy = 0, fz = 0; 
    
    // L_n^0 terms, n = 1, ..., p
    for (int n = 1; n <= p; ++n) {
      zs2 += L[midx(n, 1)] * powers_r[n - 1] * legendre[midx(n - 1, 0)] * 
        sqf[n + 1] / sqf[n - 1]; 
      fz += real(L[midx(n, 0)]) * powers_r[n - 1] * n * 
        legendre[midx(n - 1, 0)];
    }

    // L_n^m terms, n = 1, ...., p, m = 1, ..., n
    for (int n = 1; n <= p; ++n) {      
      // zs3 for z derivative
      for (int m = 1; m <= n - 1; ++m) {
        zs3 += L[midx(n, m)] * powers_ephi[m] * powers_r[n - 1] * 
          legendre[midx(n - 1, m)] * sqf[n - m] / sqf[n + m - 1] * 
          sqf[n + m] / sqf[n + m - 1];
      }

      // zs2 for x, y derivatives
      for (int m = 2; m <= n; ++m) {
        zs2 += L[midx(n, m)] * powers_ephi[m - 1] * powers_r[n - 1] * 
          legendre[midx(n - 1, m - 1)] * sqf[n - m] / sqf[n + m - 2] * 
          sqf[n + m] / sqf[n + m - 2];
      }

      // zs1 for x, y derivatives
      for (int m = 0; m <= n - 2; ++m) {
        zs1 += L[midx(n, m)] * powers_ephi[m + 1] * powers_r[n - 1] * 
          legendre[midx(n - 1, m + 1)] * sqf[n - m] / sqf[n + m];
      }
    }

    fx = real(zs2 - zs1); 
    fy = -imag(zs2 + zs1); 
    fz += 2.0 * real(zs3); 

    retval.push_back(fx * scale); 
    retval.push_back(fy * scale); 
    retval.push_back(-fz * scale);
  }

  return retval; 
}

} // namespace dashmm



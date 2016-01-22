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


/// \file src/laplace_sph.cc
/// \brief Implementation of LaplaceSPH

#include "include/laplace_sph.h"

namespace dashmm {

std::map<int, uLaplaceSPHTable> builtin_laplace_table_;

LaplaceSPH::LaplaceSPH(Point center, int n_digits) : n_digits_{n_digits} {
  LaplaceSPHTableIterator entry = builtin_laplace_table_.find(n_digits);
  assert(entry != builtin_laplace_table_.end());
  uLaplaceSPHTable &table = entry->second;
  int p = table->p();
  int n_terms = (p + 1) * (p + 2) / 2;
  bytes_ = sizeof(LaplaceSPHData) + sizeof(dcomplex_t) * n_terms;
  data_ = static_cast<LaplaceSPHData *>(malloc(bytes_));
  assert(valid());
  data_->type = type();
  data_->n_digits = n_digits;
  data_->center = center;
  for (int i = 0; i < n_terms; ++i)
    data_->expansion[i] = 0;
}

LaplaceSPH::LaplaceSPH(LaplaceSPHData *ptr, size_t bytes, int n_digits)
  : n_digits_{n_digits} {
  data_ = ptr;
  bytes_ = bytes;
  if (data_)
    data_->n_digits = n_digits;
}

LaplaceSPH::~LaplaceSPH() {
  if (valid()) {
    free(data_);
    data_ = nullptr;
  }
}

std::unique_ptr<Expansion> LaplaceSPH::S_to_M(Point center,
                                              Source *first, Source *last,
                                              double scale) const {
  LaplaceSPH *retval{new LaplaceSPH{center, n_digits_}};
  dcomplex_t *expansion = &retval->data_->expansion[0];
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();
  const double *sqf = table->sqf();

  double *legendre = new double[(p + 1) * (p + 2) / 2];
  double *powers_r = new double[p + 1];
  dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
  powers_r[0] = 1.0;
  powers_ephi[0] = dcomplex_t{1.0, 0.0};

  for (auto i = first; i != last; ++i) {
    Point dist = point_sub(i->position, center);
    double q = i->charge;
    double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
    double r = dist.norm();

    double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

    // Compute exp(-i * phi) for the azimuthal angle phi
    dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                       dcomplex_t{dist.x() / proj, -dist.y() / proj});

    // Compute powers of r
    r *= scale;
    for (int j = 1; j <= p; ++j) {
      powers_r[j] = powers_r[j - 1] * r;
    }

    // Compute powers of exp(-i * phi)
    for (int j = 1; j <= p; ++j) {
      powers_ephi[j] = powers_ephi[j - 1] * ephi;
    }

    // Compute multipole expansion M_n^m
    legendre_Plm(p, ctheta, legendre);
    for (int m = 0; m <= p; ++m) {
      for (int n = m; n <= p; ++n) {
        expansion[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
          legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
      }
    }
  }

  delete [] legendre;
  delete [] powers_r;
  delete [] powers_ephi;
  return std::unique_ptr<Expansion>{retval};
}

std::unique_ptr<Expansion> LaplaceSPH::S_to_L(Point center,
                                              Source *first, Source *last,
                                              double scale) const {
  LaplaceSPH *retval{new LaplaceSPH{center, n_digits_}};
  dcomplex_t *expansion = &retval->data_->expansion[0];
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();
  const double *sqf = table->sqf();

  double *legendre = new double[(p + 1) * (p + 2) / 2];
  double *powers_r = new double[p + 1];
  dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
  powers_ephi[0] = dcomplex_t{1.0, 0.0};

  for (auto i = first; i != last; ++i) {
    Point dist = point_sub(i->position, center);
    double q = i->charge;
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
    for (int j = 1; j <= p; ++j) {
      powers_r[j] = powers_r[j - 1] / r;
    }

    // Compute powers of exp(-i * phi)
    for (int j = 1; j <= p; ++j) {
      powers_ephi[j] = powers_ephi[j - 1] * ephi;
    }

    // compute local expansion L_n^m
    legendre_Plm(p, ctheta, legendre);
    for (int m = 0; m <= p; ++m) {
      for (int n = m; n <= p; ++n) {
        expansion[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
          legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
      }
    }
  }

  delete [] legendre;
  delete [] powers_r;
  delete [] powers_ephi;

  return std::unique_ptr<Expansion>{retval};
}

std::unique_ptr<Expansion> LaplaceSPH::M_to_M(int from_child,
                                              double s_size) const {
  // The function is called on th expansion of the child box and
  // s_size is the child box's size.
  double h = s_size / 2;
  double px = data_->center.x() + (from_child % 2 == 0 ? h : -h);
  double py = data_->center.y() + (from_child % 4 <= 1 ? h : -h);
  double pz = data_->center.z() + (from_child < 4 ? h : -h);

  LaplaceSPH *retval{new LaplaceSPH{Point{px, py, pz}, n_digits_}};

  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();
  const double *sqbinom = table->sqbinom();

  // Get precomputed Wigner d-matrix for rotation about the y-axis
  const double *d1 = (from_child < 4 ?
                      table->dmat_plus(1.0 / sqrt(3.0)) :
                      table->dmat_plus(-1.0 / sqrt(3.0)));
  const double *d2 = (from_child < 4 ?
                      table->dmat_minus(1.0 / sqrt(3.0)) :
                      table->dmat_minus(-1.0 / sqrt(3.0)));

  // Shift distance along the z-axis, combined with Y_n^0(pi, 0)
  const double rho = -sqrt(3) / 2;

  // Compute powers of rho
  double *powers_rho = new double[p + 1];
  powers_rho[0] = 1.0;
  for (int i = 1; i <= p; ++i) {
    powers_rho[i] = powers_rho[i - 1] * rho;
  }

  dcomplex_t *W1 = &retval->data_->expansion[0];
  dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

  // Table of rotation angle about the z-axis, as an integer multiple of pi / 4
  const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
  // Get rotation angle about the z-axis
  double alpha = tab_alpha[from_child] * M_PI_4;

  // Rotate the multipole expansion of the child box about z-axis
  dcomplex_t *M = &data_->expansion[0];
  retval->rotate_sph_z(M, alpha, W1);

  // Rotate the previous result further about the y-axis
  retval->rotate_sph_y(W1, d1, W2);

  // Offset to operate multipole expansion
  int offset = 0;

  // Shift along the z-axis by a distance of rho, write result in W1
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      W1[offset] = W2[offset];
      for (int k = 1; k <= n - m; ++k) {
        W1[offset] += W2[midx(n - k, m)] * powers_rho[k] *
        sqbinom[midx(n - m, k)] * sqbinom[midx(n + m, k)];
      }
      offset++;
    }
  }

  // Reverse rotate the shifted harmonic expansion about the y-axis
  retval->rotate_sph_y(W1, d2, W2);

  // Reverse rotate the previous result further about the z-axis
  retval->rotate_sph_z(W2, -alpha, W1);

  double scale = 1;
  offset = 0;
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      W1[offset++] *= scale;
    }
    scale /= 2;
  }

  delete [] W2;
  delete [] powers_rho;
  return std::unique_ptr<Expansion>{retval};
}

std::unique_ptr<Expansion> LaplaceSPH::M_to_L(Index s_index, double s_size,
                                              Index t_index) const {
  int t2s_x = s_index.x() - t_index.x();
  int t2s_y = s_index.y() - t_index.y();
  int t2s_z = s_index.z() - t_index.z();
  double tx = data_->center.x() - t2s_x * s_size;
  double ty = data_->center.y() - t2s_y * s_size;
  double tz = data_->center.z() - t2s_z * s_size;

  LaplaceSPH *retval{new LaplaceSPH{Point{tx, ty, tz}, n_digits_}};
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();


  // Shifting distance
  double rho = sqrt(t2s_x * t2s_x + t2s_y * t2s_y + t2s_z * t2s_z);

  // Compute powers of rho
  double *powers_rho = new double[p * 2 + 1];
  powers_rho[0] = 1.0 / rho;
  for (int i = 1; i <= p * 2; i++) {
    powers_rho[i] = powers_rho[i - 1] / rho;
  }

  // Scaling factor
  double scale = 1.0 / s_size;

  // Temporary space to hold rotated spherical harmonic
  dcomplex_t *W1 = &retval->data_->expansion[0];
  dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

  // Compute the projection of t2s on the x-y plane
  const double proj = sqrt(t2s_x * t2s_x + t2s_y * t2s_y);

  // Handle of the multipole expansion
  dcomplex_t *M = &data_->expansion[0];

  if (proj < 1e-14) {
    if (t2s_z > 0) {
      retval->M_to_L_zp(M, powers_rho, scale, W1);
    } else {
      retval->M_to_L_zm(M, powers_rho, scale, W1);
    }
  } else {
    // azimuthal angle
    double beta = acos(t2s_x / proj);
    if (t2s_y < 0) {
      beta = 2 * M_PI - beta;
    }

    // Get precomputed Wigner d-matrix for rotation about y-axis
    const double *d1 = table->dmat_plus(t2s_z / rho);
    const double *d2 = table->dmat_minus(t2s_z / rho);
    retval->rotate_sph_z(M, beta, W1);
    retval->rotate_sph_y(W1, d1, W2);
    retval->M_to_L_zp(W2, powers_rho, scale, W1);
    retval->rotate_sph_y(W1, d2, W2);
    retval->rotate_sph_z(W2, -beta, W1);
  }

  delete [] W2;
  delete [] powers_rho;
  return std::unique_ptr<Expansion>{retval};
}

std::unique_ptr<Expansion> LaplaceSPH::L_to_L(int to_child,
                                              double t_size) const {
  // The function is called on the parent box and t_size is its child size
  double h = t_size / 2;
  double cx = data_->center.x() + (to_child % 2 == 0 ? -h : h);
  double cy = data_->center.y() + (to_child % 4 <= 1 ? -h : h);
  double cz = data_->center.z() + (to_child < 4 ? -h : h);

  LaplaceSPH *retval{new LaplaceSPH{Point{cx, cy, cz}, n_digits_}};

  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();
  const double *sqbinom = table->sqbinom();

  // Get precomputed Wigner d-matrix for rotation about the y-axis
  const double *d1 = (to_child < 4 ?
                      table->dmat_plus(1.0 / sqrt(3.0)) :
                      table->dmat_plus(-1.0 / sqrt(3.0)));
  const double *d2 = (to_child < 4 ?
                      table->dmat_minus(1.0 / sqrt(3.0)) :
                      table->dmat_minus(-1.0 / sqrt(3.0)));

  // Shift distance along the z-axis, combined with Y_n^0(pi, 0)
  const double rho = -sqrt(3) / 4;

  // Compute powers of rho
  double *powers_rho = new double[p + 1];
  powers_rho[0] = 1.0;
  for (int i = 1; i <= p; ++i) {
    powers_rho[i] = powers_rho[i - 1] * rho;
  }

  dcomplex_t *W1 = &retval->data_->expansion[0];
  dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

  // Table of rotation angle about the z-axis as an integer multiple of pi / 4
  const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
  // Get rotation angle about the z-axis
  double alpha = tab_alpha[to_child] * M_PI_4;

  // Rotate the local expansion of the parent box about z-axis
  dcomplex_t *L = &data_->expansion[0];
  retval->rotate_sph_z(L, alpha, W1);

  // Rotate the previous result further about the y-axis
  retval->rotate_sph_y(W1, d1, W2);

  // Offset to operate local expansion
  int offset = 0;

  // Shift along the z-axis by a distance of rho, write result in W1
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      W1[offset] = W2[offset];
      for (int k = 1; k <= p - n; k++) {
        W1[offset] += W2[midx(n + k, m)] * powers_rho[k] *
          sqbinom[midx(n + k - m, k)] * sqbinom[midx(n + k + m, k)];
      }
      offset++;
    }
  }

  // Reverse rotate the shifted harmonic expansion about the y-axis
  retval->rotate_sph_y(W1, d2, W2);

  // Reverse rotate the previous result further about the z-axis
  retval->rotate_sph_z(W2, -alpha, W1);

  double scale = 1;
  offset = 0;
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      W1[offset++] *= scale;
    }
    scale /= 2;
  }

  delete [] W2;
  delete [] powers_rho;
  return std::unique_ptr<Expansion>{retval};
}

void LaplaceSPH::M_to_T(Target *first, Target *last, double scale) const {
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();
  const double *sqf = table->sqf();

  double *legendre = new double[(p + 1) * (p + 2) / 2];
  double *powers_r = new double[p + 1];
  dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
  powers_ephi[0] = dcomplex_t{1.0, 0.0};
  dcomplex_t *M = &data_->expansion[0];

  for (auto i = first; i != last; ++i) {
    Point dist = point_sub(i->position, data_->center);
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
    for (int j = 1; j <= p; ++j) {
      powers_r[j] = powers_r[j - 1] / r;
    }

    // Compute powers of exp(i * phi)
    for (int j = 1; j <= p; ++j) {
      powers_ephi[j] = powers_ephi[j - 1] * ephi;
    }

    // Evaluate the multipole expansion M_n^0
    legendre_Plm(p, ctheta, legendre);
    for (int n = 0; n <= p; ++n) {
      potential += M[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)];
    }

    // Evaluate the multipole expansions M_n^m, where m = 1, ..., p
    for (int m = 1; m <= p; ++m) {
      for (int n = m; n <= p; ++n) {
        potential += 2.0 * real(M[midx(n, m)] * powers_ephi[m]) *
          powers_r[n] * legendre[midx(n, m)] * sqf[n - m] /sqf[n + m];
      }
    }

    i->phi += potential;
  }

  delete [] powers_r;
  delete [] powers_ephi;
  delete [] legendre;
}

void LaplaceSPH::L_to_T(Target *first, Target *last, double scale) const {
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();
  const double *sqf = table->sqf();

  double *legendre = new double[(p + 1) * (p + 2) / 2];
  double *powers_r = new double[p + 1];
  dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
  powers_r[0] = 1.0;
  powers_ephi[0] = dcomplex_t{1.0, 0.0};

  dcomplex_t *L = &data_->expansion[0];

  for (auto i = first; i != last; ++i) {
    Point dist = point_sub(i->position, data_->center);
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
    for (int j = 1; j <= p; ++j) {
      powers_r[j] = powers_r[j - 1] * r;
    }

    // Compute powers of exp(i * phi)
    for (int j = 1; j <= p; ++j) {
      powers_ephi[j] = powers_ephi[j - 1] * ephi;
    }

    // Evaluate the local expansion L_n^0
    legendre_Plm(p, ctheta, legendre);
    for (int n = 0; n <= p; ++n) {
      potential += L[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)];
    }

    // Evaluate the local expansions L_n^m, where m = 1, ..., p
    for (int m = 1; m <= p; ++m) {
      for (int n = m; n <= p; ++n) {
        potential += 2.0 * real(L[midx(n, m)] * powers_ephi[m]) *
          powers_r[n] * legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
      }
    }

    i->phi += potential;
  }

  delete [] powers_r;
  delete [] powers_ephi;
  delete [] legendre;
}

void LaplaceSPH::S_to_T(Source *s_first, Source *s_last,
                        Target *t_first, Target *t_last) const {
  for (auto i = t_first; i != t_last; ++i) {
    dcomplex_t potential{0.0, 0.0};
    for (auto j = s_first; j != s_last; ++j) {
      Point s2t = point_sub(i->position, j->position);
      double dist = s2t.norm();
      if (dist > 0)
        potential += j->charge / dist;
    }
    i->phi += potential;
  }
}

void LaplaceSPH::add_expansion(const Expansion *temp1) {
  dcomplex_t *expansion = &data_->expansion[0];
  for (size_t i = 0; i < temp1->size(); ++i)
    expansion[i] += temp1->term(i);
}

std::unique_ptr<Expansion> LaplaceSPH::get_new_expansion(Point center) const {
  LaplaceSPH *retval{new LaplaceSPH{center, n_digits_}};
  return std::unique_ptr<Expansion>{retval};
}

void LaplaceSPH::rotate_sph_z(const dcomplex_t *M, double alpha,
                              dcomplex_t *MR) {
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();

  // Compute exp(i * alpha)
  dcomplex_t ealpha{cos(alpha), sin(alpha)};

  // Compute powers of exp(i * alpha)
  dcomplex_t *powers_ealpha = new dcomplex_t[p + 1];
  powers_ealpha[0] = dcomplex_t{1.0, 0.0};
  for (int j = 1; j <= p; ++j) {
    powers_ealpha[j] = powers_ealpha[j - 1] * ealpha;
  }

  int offset = 0;
  for (int n = 0; n <= p; n++) {
    for (int m = 0; m <= n; m++) {
      MR[offset] = M[offset] * powers_ealpha[m];
      offset++;
    }
  }

  delete [] powers_ealpha;
}

void LaplaceSPH::rotate_sph_y(const dcomplex_t *M, const double *d,
                              dcomplex_t *MR) {
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();

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

void LaplaceSPH::M_to_L_zp(const dcomplex_t *M, const double *rho,
                           double scale, dcomplex_t *L) {
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();
  const double *sqbinom = table->sqbinom();

  int offset = 0;
  for (int j = 0; j <= p; ++j) {
    for (int k = 0; k <= j; ++k) {
      L[offset] = 0;
      for (int n = k; n <= p; ++n) {
        L[offset] += M[midx(n, k)] * pow_m1(n + k) * rho[j + n] *
          sqbinom[midx(n + j, n - k)] * sqbinom[midx(n + j, n + k)];
      }
      L[offset] *= scale;
      offset++;
    }
  }
}

void LaplaceSPH::M_to_L_zm(const dcomplex_t *M, const double *rho,
                           double scale, dcomplex_t *L) {
  uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
  int p = table->p();
  const double *sqbinom = table->sqbinom();

  int offset = 0;
  for (int j = 0; j <= p; ++j) {
    for (int k = 0; k <= j; ++k) {
      L[offset] = 0;
      for (int n = k; n <= p; ++n) {
        L[offset] += M[midx(n, k)] * pow_m1(k + j) * rho[j + n] *
          sqbinom[midx(n + j, n - k)] * sqbinom[midx(n + j, n + k)];
      }
      L[offset] *= scale;
      offset++;
    }
  }
}

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

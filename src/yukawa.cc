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
/// \brief Implementation of Yukawa kernel


#include "builtins/yukawa.h"


namespace dashmm {

void yuk_rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR) {
  int p = builtin_yukawa_table_->p();
  // Compute exp(i * alpha)
  dcomplex_t ealpha = dcomplex_t{cos(alpha),  sin(alpha)};

  // Compute powers of exp(i * alpha)
  std::vector<dcomplex_t> powers_ealpha(p + 1); 
  powers_ealpha[0] = dcomplex_t{1.0, 0.0};
  for (int j = 1; j <= p; ++j) 
    powers_ealpha[j] = powers_ealpha[j - 1] * ealpha;
  
  int offset = 0;
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      MR[offset] = M[offset] * powers_ealpha[m];
      offset++;
    }
  }
}

void yuk_rotate_sph_y(const dcomplex_t *M, const double *d, dcomplex_t *MR) {
  int p = builtin_yukawa_table_->p();

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

void yuk_s_to_m(Point dist, double q, double scale, dcomplex_t *M) {
  int p = builtin_yukawa_table_->p();
  const double *sqf = builtin_yukawa_table_->sqf();
  double lambda = builtin_yukawa_table_->lambda();

  std::vector<double> legendre((p + 1) * (p + 2) / 2); 
  std::vector<dcomplex_t> powers_ephi(p + 1); 
  std::vector<double> bessel(p + 1); 

  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();
  
  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);

  // Compute exp(-i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, -dist.y() /proj});

  // Compute powers of exp(-i * phi)
  powers_ephi[0] = 1.0;
  for (int j = 1; j <= p; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;
  
  // Compute scaled modified spherical bessel function
  bessel_in_scaled(p, lambda * r, scale, bessel.data());

  // Compute legendre polynomial
  legendre_Plm(p, ctheta, legendre.data());

  // Compute multipole expansion M_n^m
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      int idx = midx(n, m);
      M[idx] += q * bessel[n] * sqf[idx] * legendre[idx] * powers_ephi[m];
    }
  }
}

void yuk_s_to_l(Point dist, double q, double scale, dcomplex_t *L) {
  int p = builtin_yukawa_table_->p();
  const double *sqf = builtin_yukawa_table_->sqf();
  double lambda = builtin_yukawa_table_->lambda();

  std::vector<double> legendre((p + 1) * (p + 2) / 2); 
  std::vector<double> bessel(p + 1); 
  std::vector<dcomplex_t> powers_ephi(p + 1); 

  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();

  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);
  
  // Compute exp(-i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, -dist.y() / proj});
  
  // Compute powers of exp(-i * phi)
  powers_ephi[0] = 1.0;
  for (int j = 1; j <= p; j++) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;
      
  // Compute scaled modified spherical bessel function
  bessel_kn_scaled(p, lambda * r, scale, bessel.data());

  // Compute legendre polynomial
  legendre_Plm(p, ctheta, legendre.data());

  // Compute local expansion L_n^m
  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      int idx = midx(n, m);
      L[idx] += q * bessel[n] * sqf[idx] * legendre[idx] * powers_ephi[m];
    }
  }
}

void yuk_m_to_m(int from_child, const dcomplex_t *M, double scale, 
                dcomplex_t *W) {
  int p = builtin_yukawa_table_->p();

  // Get precomputed Wigner d-matrix for rotation about the y-axis
  const double *d1 = (from_child < 4 ?
                      builtin_yukawa_table_->dmat_plus(1.0 / sqrt(3.0)) :
                      builtin_yukawa_table_->dmat_plus(-1.0 / sqrt(3.0)));
  const double *d2 = (from_child < 4 ?
                      builtin_yukawa_table_->dmat_minus(1.0 / sqrt(3.0)) :
                      builtin_yukawa_table_->dmat_minus(-1.0 / sqrt(3.0)));
  
  // Get precomputed coefficients for shifting along z-axis
  const double *coeff = builtin_yukawa_table_->m2m(scale);

  // Table of rotation angle about z-axis, as an integer of multiple of pi / 4
  const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};

  // Get rotation angle
  double alpha = tab_alpha[from_child] * M_PI_4;

  std::vector<dcomplex_t> T((p + 1) * (p + 2) / 2); 

  yuk_rotate_sph_z(M, alpha, W);
  yuk_rotate_sph_y(W, d1, T.data());

  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      dcomplex_t temp{0.0, 0.0};
      for (int np = m; np <= p; ++np) {
        temp += T[midx(np, m)] * coeff[sidx(n, m, np, p)];
      }
      W[midx(n, m)] = temp;
    }
  }

  yuk_rotate_sph_y(W, d2, T.data());
  yuk_rotate_sph_z(T.data(), -alpha, W);
}

void yuk_l_to_l(int to_child, const dcomplex_t *L, double scale, 
                dcomplex_t *W) { 
  int p = builtin_yukawa_table_->p();

  // Get precomputed Wigner d-matrix for rotation about the y-axis
  const double *d1 = (to_child < 4 ?
                      builtin_yukawa_table_->dmat_plus(1.0 / sqrt(3)) :
                      builtin_yukawa_table_->dmat_plus(-1.0 / sqrt(3)));
  const double *d2 = (to_child < 4 ?
                      builtin_yukawa_table_->dmat_minus(1.0 / sqrt(3)) :
                      builtin_yukawa_table_->dmat_minus(-1.0 / sqrt(3)));
  
  // Get precomputed coefficients for shifting along z-axis
  const double *coeff = builtin_yukawa_table_->l2l(scale);

  // Table of rotation angle about z-axis, as an integer multiple of pi / 4
  const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
  
  // Get rotation angle
  double alpha = tab_alpha[to_child] * M_PI_4;

  std::vector<dcomplex_t> T((p + 1) * (p + 2) / 2); 

  yuk_rotate_sph_z(L, alpha, W);
  yuk_rotate_sph_y(W, d1, T.data());

  for (int n = 0; n <= p; ++n) {
    for (int m = 0; m <= n; ++m) {
      dcomplex_t temp{0.0, 0.0};
      for (int np = m; np <= p; ++np) {
        temp += T[midx(np, m)] * coeff[sidx(n, m, np, p)];
      }
      W[midx(n, m)] = temp;
    }
  }

  yuk_rotate_sph_y(W, d2, T.data());
  yuk_rotate_sph_z(T.data(), -alpha, W);
}


std::vector<double> yuk_m_to_t(Point dist, double scale, 
                               const dcomplex_t *M, bool g) {
  std::vector<double> retval; 

  int p = builtin_yukawa_table_->p();
  double lambda = builtin_yukawa_table_->lambda();
  std::vector<double> legendre((p + 2) * (p + 3) / 2); 
  std::vector<dcomplex_t> temp((p + 2) * (p + 3) / 2); 
  std::vector<double> bessel(p + 2);  
  std::vector<dcomplex_t> powers_ephi(p + 2); 

  // Compute potential first 
  dcomplex_t potential{0.0, 0.0};
  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();

  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);
  
  // Compute exp(i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, dist.y() / proj});
  
  // Compute powers of exp(i * phi)
  powers_ephi[0] = 1.0;
  for (int j = 1; j <= p + 1; j++) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;

  // Compute scaled modified spherical bessel function
  bessel_kn_scaled(p + 1, lambda * r, scale, bessel.data());

  // Compute legendre polynomial
  legendre_Plm(p + 1, ctheta, legendre.data());
  
  // Evaluate M_n^0
  for (int n = 0; n <= p; ++n) {
    temp[midx(n, 0)] = bessel[n] * legendre[midx(n, 0)]; 
    potential += M[midx(n, 0)] * temp[midx(n, 0)]; 
 }
  
  // Evaluate M_n^m
  for (int n = 1; n <= p; ++n) {
    for (int m = 1; m <= n; ++m) {
      temp[midx(n, m)] = bessel[n] * legendre[midx(n, m)] * powers_ephi[m]; 
      potential += 2.0 * real(M[midx(n, m)] * temp[midx(n, m)]);
    }
  }

  retval.push_back(real(potential)); 

  if (g) {
    for (int m = 0; m <= p + 1; ++m) {
      temp[midx(p + 1, m)] = bessel[p + 1] * legendre[midx(p + 1, m)] * 
        powers_ephi[m];
    }

    // Compute field values 
    dcomplex_t rpotz = 0.0, cpz = 0.0, tx = 0, ty = 0, tz = 0; 
    double fx = 0, fy = 0, fz = 0;

    // Contribution from n = 0, m = 0 term 
    rpotz = M[midx(0, 0)] / scale; 
    cpz = rpotz * temp[midx(1, 1)]; 
    tx = -conj(cpz); 
    ty = -rpotz * temp[midx(1, 0)]; 
    tz = cpz; 

    // Contribution from n = 1, m = 0 term 
    rpotz = M[midx(1, 0)] / 3.0; 
    cpz = rpotz / scale * temp[midx(2, 1)]; 
    tx -= conj(cpz); 
    ty -= rpotz * (temp[midx(0, 0)] * scale + 2.0 * temp[midx(2, 0)] / scale);
    tz += cpz;

    // Contribution from n > 1, m = 0 terms 
    for (int n = 2; n <= p; ++n) {
      rpotz = M[midx(n, 0)] / (2.0 * n + 1); 
      cpz = rpotz * (temp[midx(n - 1, 1)] * scale - 
                     temp[midx(n + 1, 1)] / scale);
      tx += conj(cpz);
      ty -= rpotz * (real(temp[midx(n - 1, 0)]) * ((double) n) * scale + 
                     real(temp[midx(n + 1, 0)]) * ((double) (n + 1)) / scale);
      tz -= cpz;
    }
    
    fx = -lambda * real(tx - tz) / 2.0; 
    fy = -lambda * real((tx + tz) * dcomplex_t{0.0, 1.0}) / 2.0; 
    fz += lambda * real(ty);

    // Contribution from n = 1, m = 1
    cpz = M[midx(1, 1)] / 3.0; 
    tx = -cpz * (temp[midx(0, 0)] * scale - temp[midx(2, 0)] / scale) * 2.0; 
    ty = -cpz * temp[midx(2, 1)] / scale; 
    tz = cpz * temp[midx(2, 2)] / scale; 

    // Contributions for n = 2, ..., p
    for (int n = 2; n <= p; ++n) {
      for (int m = 1; m <= n - 2; ++m) {
        cpz = M[midx(n, m)] / (2.0 * n + 1); 
        tx -= cpz * (temp[midx(n - 1, m - 1)] * 
                     ((double) (n + m - 1) * (n + m)) * scale - 
                     temp[midx(n + 1, m - 1)] * 
                     ((double) (n - m + 1) * (n - m + 2)) / scale);
        ty -= cpz * (temp[midx(n - 1, m)] * ((double) (n + m)) * scale + 
                     temp[midx(n + 1, m)] * ((double) (n - m + 1)) / scale);
        tz -= cpz * (temp[midx(n - 1, m + 1)] * scale - 
                     temp[midx(n + 1, m + 1)] / scale); 
      }

      // m = n - 1
      cpz = M[midx(n, n - 1)] / (2.0 * n + 1); 
      tx -= cpz * (temp[midx(n - 1, n - 2)] * 
                   ((double) (n + n - 2) * (n + n - 1)) * scale - 
                   temp[midx(n + 1, n - 2)] * 6.0 / scale);
      ty -= cpz * (temp[midx(n - 1, n - 1)] * ((double) (n + n - 1)) * scale + 
                   temp[midx(n + 1, n - 1)] * 2.0 / scale );
      tz += cpz * temp[midx(n + 1, n)] / scale;

      // m = n
      cpz = M[midx(n, n)] / (2.0 * n + 1); 
      tx -= cpz * (temp[midx(n - 1, n - 1)] * 
                   ((double) (n + n - 1) * (n + n)) * scale - 
                   temp[midx(n + 1, n - 1)] * 2.0 / scale);
      ty -= cpz * temp[midx(n + 1, n)] / scale; 
      tz += cpz * temp[midx(n + 1, n + 1)] / scale; 
    }

    fx -= lambda * real(tx - tz); 
    fy -= lambda * real((tx + tz) * dcomplex_t{0.0, 1.0}); 
    fz += 2.0 * lambda * real(ty);

    retval.push_back(fx); 
    retval.push_back(fy);
    retval.push_back(fz); 
  }
  return retval;
}

std::vector<double> yuk_l_to_t(Point dist, double scale, 
                               const dcomplex_t *L, bool g) {
  std::vector<double> retval; 

  int p = builtin_yukawa_table_->p();
  double lambda = builtin_yukawa_table_->lambda();
  std::vector<double> legendre((p + 2) * (p + 3) / 2);
  std::vector<dcomplex_t> temp((p + 2) * (p + 3) / 2); 
  std::vector<double> bessel(p + 2); 
  std::vector<dcomplex_t> powers_ephi(p + 2); 

  // Compute potential first
  dcomplex_t potential{0.0, 0.0};
  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();
  
  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);
  
  // Compute exp(i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, dist.y() / proj});
  
  // Compute powers of exp(i * phi)
  powers_ephi[0] = 1.0;
  for (int j = 1; j <= p + 1; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;

  // Compute scaled modified spherical bessel function
  bessel_in_scaled(p + 1, lambda * r, scale, bessel.data());

  // Compute legendre polynomial
  legendre_Plm(p + 1, ctheta, legendre.data());

  // Evaluate local expansion L_n^0
  for (int n = 0; n <= p; ++n) {
    temp[midx(n, 0)] = bessel[n] * legendre[midx(n, 0)]; 
    potential += L[midx(n, 0)] * temp[midx(n, 0)]; 
  }
      
  // Evaluate L_n^m
  for (int n = 1; n <= p; ++n) {
    for (int m = 1; m <= n; ++m) {
      temp[midx(n, m)] = bessel[n] * legendre[midx(n, m)] * powers_ephi[m];
      potential += 2.0 * real(L[midx(n, m)] * temp[midx(n, m)]); 
    }
  }

  retval.push_back(real(potential)); 

  if (g) {
    for (int m = 0; m <= p + 1; ++m) {
      temp[midx(p + 1, m)] = bessel[p + 1] * legendre[midx(p + 1, m)] * 
        powers_ephi[m];
    }

    // Compute field values
    dcomplex_t rpotz = 0.0, cpz = 0, tx = 0, ty = 0, tz = 0; 
    double fx = 0, fy = 0, fz = 0;

    // Contribution from n = 0, m = 0 term 
    rpotz = L[midx(0, 0)] * scale; 
    cpz = rpotz * temp[midx(1, 1)]; 
    tx = conj(cpz); 
    ty = rpotz * temp[midx(1, 0)]; 
    tz = -cpz; 

    // Contribution from n = 1, m = 0 term 
    rpotz = L[midx(1, 0)] / 3.0; 
    cpz = rpotz * scale * temp[midx(2, 1)]; 
    tx += conj(cpz); 
    ty += rpotz * (temp[midx(0, 0)] / scale + 
                   2.0 * temp[midx(2, 0)] * scale); 
    tz -= cpz;
    
    // Contribution from n > 1, m = 0 term 
    for (int n = 2; n <= p; ++n) {
      rpotz = L[midx(n, 0)] / (2.0 * n + 1); 
      cpz = rpotz * (temp[midx(n - 1, 1)] / scale - 
                     temp[midx(n + 1, 1)] * scale);
      tx -= conj(cpz); 
      ty += rpotz * (real(temp[midx(n - 1, 0)]) * ((double) n) / scale + 
                     real(temp[midx(n + 1, 0)]) * ((double) (n + 1)) * scale);
      tz += cpz;
    }

    fx -= lambda * real(tx - tz) / 2.0; 
    fy -= lambda * real((tx + tz) * dcomplex_t{0.0, 1.0}) / 2.0; 
    fz += lambda * real(ty); 

    // Contribution from n = 1, m = 1
    cpz = L[midx(1, 1)] / 3.0; 
    tx = cpz * (real(temp[midx(0, 0)]) / scale - 
                real(temp[midx(2, 0)]) * scale) * 2.0; 
    cpz *= scale; 
    ty = cpz * temp[midx(2, 1)]; 
    tz = -cpz * temp[midx(2, 2)]; 
    
    // Contribution from n = 2, ..., p, 
    for (int n = 2; n <= p; ++n) {
      for (int m = 1; m <= n - 2; ++m) {
        cpz = L[midx(n, m)] / (2.0 * n + 1); 
        tx += cpz * (temp[midx(n - 1, m - 1)] * 
                     ((double) (n + m - 1) * (n + m)) / scale -
                     temp[midx(n + 1, m - 1)] * 
                     ((double) (n - m + 1) * (n - m + 2)) * scale);
        ty += cpz * (temp[midx(n - 1, m)] * ((double) (n + m)) / scale + 
                     temp[midx(n + 1, m)] * ((double) (n - m + 1)) * scale);
        tz += cpz * (temp[midx(n - 1, m + 1)] / scale - 
                     temp[midx(n + 1, m + 1)] * scale);
      }
      
      // m = n - 1
      cpz = L[midx(n, n - 1)] / (2.0 * n + 1); 
      tx += cpz * (temp[midx(n - 1, n - 2)] * 
                   ((double) (n + n - 2) * (n + n - 1)) / scale - 
                   temp[midx(n + 1, n - 2)] * 6.0 * scale);
      ty += cpz * (temp[midx(n - 1, n - 1)] * 
                   ((double) (n + n - 1)) / scale + 
                   temp[midx(n + 1, n - 1)] * 2.0 * scale);
      tz -= cpz * temp[midx(n + 1, n)] * scale;

      // m = n
      cpz = L[midx(n, n)] / (2.0 * n + 1); 
      tx += cpz * (temp[midx(n - 1, n - 1)] * 
                   ((double) (n + n - 1) * (n + n)) / scale - 
                   temp[midx(n + 1, n - 1)] * 2.0 * scale); 
      ty += cpz * temp[midx(n + 1, n)] * scale; 
      tz -= cpz * temp[midx(n + 1, n + 1)] * scale; 
    }

    fx -= lambda * real(tx - tz); 
    fy -= lambda * real((tx + tz) * dcomplex_t{0.0, 1.0}); 
    fz += 2.0 * lambda * real(ty);

    retval.push_back(fx); 
    retval.push_back(fy);
    retval.push_back(fz); 
  }

  return retval; 
}

void yuk_m_to_i(const dcomplex_t *M, ViewSet &views, double scale, int id) {
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
  
  int p = builtin_yukawa_table_->p();
  double ld = builtin_yukawa_table_->lambda() * 
    builtin_yukawa_table_->size(scale);

  // Get precomputed Wigner d-matrix
  const double *d1 = builtin_yukawa_table_->dmat_plus(0.0);
  const double *d2 = builtin_yukawa_table_->dmat_minus(0.0);
  
  // Get number of Gaussian quadrature points
  int s = builtin_yukawa_table_->s();
  
  // Get Gaussian quadrature points
  const double *x = builtin_yukawa_table_->x();
  
  // Get number of Fourier modes
  const int *f = builtin_yukawa_table_->f();
  
  // Get number of trapezoid points
  const int *mk = builtin_yukawa_table_->m(scale);
  const int *smf = builtin_yukawa_table_->smf(scale);
  
  // Get e^{i * m * alpha_j}
  const dcomplex_t *ealphaj = builtin_yukawa_table_->ealphaj(scale);
  
  // Allocate temporary space to handle x-/y-direction expansion
  std::vector<dcomplex_t> W1((p + 1) * (p + 2) / 2); 
  std::vector<dcomplex_t> W2((p + 1) * (p + 2) / 2); 

  // Setup y-direction
  yuk_rotate_sph_z(M, -M_PI / 2, W1.data());
  yuk_rotate_sph_y(W1.data(), d2, W2.data());

  // Setup x-direction
  yuk_rotate_sph_y(M, d1, W1.data());

  // Addresses of the spherical harmonic expansions
  const dcomplex_t *SH[3] = {W1.data(), W2.data(), M};

  std::vector<double> legendre((p + 1) * (p + 2) / 2);

  for (int dir = 0; dir <=2; ++dir) {
    int offset = 0;
    for (int k = 0; k < s; ++k) {
      legendre_Plm_gt1_scaled(p, 1 + x[k] / ld, scale, legendre.data());

      // Handle M_n^m where n is even
      std::vector<dcomplex_t> z1(f[k] + 1); 
      // Handle M_n^m where n is odd
      std::vector<dcomplex_t> z2(f[k] + 1); 

      // Process M_n^0 terms
      z1[0] = 0;
      z2[0] = 0;
      for (int n = 0; n <= p; n += 2) {
        z1[0] += SH[dir][midx(n, 0)] * legendre[midx(n, 0)];
      }
      for (int n = 1; n <= p; n += 2) {
        z2[0] += SH[dir][midx(n, 0)] * legendre[midx(n, 0)];
      }
      
      // Process M_n^m, where m is odd
      for (int m = 1; m <= f[k]; m += 2) {
        z1[m] = 0;
        z2[m] = 0;
        for (int n = m; n <= p; n += 2) {
          z2[m] += SH[dir][midx(n, m)] * legendre[midx(n, m)];
        }
        for (int n = m + 1; n <= p; n += 2) {
          z1[m] += SH[dir][midx(n, m)] * legendre[midx(n, m)];
        }
      }

      // Process M_n^m, where m is even
      for (int m = 2; m <= f[k]; m += 2) {
        z1[m] = 0;
        z2[m] = 0;
        for (int n = m; n <= p; n += 2) {
          z1[m] += SH[dir][midx(n, m)] * legendre[midx(n, m)];
        }
        for (int n = m + 1; n <= p; n += 2) {
          z2[m] += SH[dir][midx(n, m)] * legendre[midx(n, m)];
        }
      }
      
      // Compute W(k, j)
      for (int j = 1; j <= mk[k] / 2; ++j) {
        dcomplex_t up{z1[0] + z2[0]};
        dcomplex_t dn{z1[0] - z2[0]};
        dcomplex_t power_I{0.0, 1.0};
        for (int m = 1; m <= f[k]; ++m) {
          int idx = smf[k] + (j - 1) * f[k] + m - 1;
          up += 2 * real(ealphaj[idx] * (z1[m] + z2[m])) * power_I;
          dn += 2 * real(ealphaj[idx] * (z1[m] - z2[m])) * power_I;
          power_I *= dcomplex_t{0.0, 1.0};
        }
        EP[dir][offset] = up;
        EM[dir][offset] = dn;
        offset++;
      }      
    }
  }
}

void yuk_i_to_i(Index s_index, Index t_index, const ViewSet &s_views, 
                int sid, int tid, double scale, ViewSet &t_views) {
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
  int nexp = builtin_yukawa_table_->nexp(scale);

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

    if (tag == -1) {
      break;
    }
    
    if (tag <= 1) {
      yuk_e_to_e(T[i], S[5], dx, dy, 0, scale);
    } else if (tag <= 5) {
      yuk_e_to_e(T[i], S[3], dz, dx, 0, scale);
    } else if (tag <= 13) {
      yuk_e_to_e(T[i], S[1], -dz, dy, 0, scale);
    } else if (tag <= 15) {
      yuk_e_to_e(T[i], S[4], -dx, -dy, 0, scale);
    } else if (tag <= 19) {
      yuk_e_to_e(T[i], S[2], -dz, -dx, 0, scale);
    } else {
      yuk_e_to_e(T[i], S[0], dz, -dy, 0, scale);
    }
    
    t_views.add_view(tid + tag, view_size, C[i]);
    used[i] = true;
  }

  if (used[1] == false) {
    delete [] T2;
  }

  if (used[2] == false) {
    delete [] T3;
  }
}

void yuk_i_to_l(const ViewSet &views, int id, Index t_index, double scale, 
                dcomplex_t *L) {
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


  // t_index and t_size is the index and size of the child
  int to_child = 4 * (t_index.z() % 2) + 2 * (t_index.y() % 2) +
    (t_index.x() % 2);

  int nexp = builtin_yukawa_table_->nexp(scale);

  dcomplex_t *S = new dcomplex_t[nexp * 6]();
  dcomplex_t *S_mz = S;
  dcomplex_t *S_pz = S + nexp;
  dcomplex_t *S_my = S + 2 * nexp;
  dcomplex_t *S_py = S + 3 * nexp;
  dcomplex_t *S_mx = S + 4 * nexp;
  dcomplex_t *S_px = S + 5 * nexp;

  switch (to_child) {
  case 0:
    yuk_e_to_e(S_mz, E[uall], 0, 0, 3, scale);
    yuk_e_to_e(S_mz, E[u1234], 0, 0, 2, scale);
    yuk_e_to_e(S_pz, E[dall], 0, 0, 2, scale);
    
    yuk_e_to_e(S_my, E[nall], 0, 0, 3, scale);
    yuk_e_to_e(S_my, E[n1256], 0, 0, 2, scale);
    yuk_e_to_e(S_my, E[n12], 0, 0, 2, scale);
    yuk_e_to_e(S_py, E[sall], 0, 0, 2, scale);
    
    yuk_e_to_e(S_mx, E[eall], 0, 0, 3, scale);
    yuk_e_to_e(S_mx, E[e1357], 0, 0, 2, scale);
    yuk_e_to_e(S_mx, E[e13], 0, 0, 2, scale);
    yuk_e_to_e(S_mx, E[e1], 0, 0, 2, scale);
    yuk_e_to_e(S_px, E[wall], 0, 0, 2, scale);
    break;
  case 1:
    yuk_e_to_e(S_mz, E[uall], -1, 0, 3, scale);
    yuk_e_to_e(S_mz, E[u1234], -1, 0, 2, scale);
    yuk_e_to_e(S_pz, E[dall], 1, 0, 2, scale);
    
    yuk_e_to_e(S_my, E[nall], 0, -1, 3, scale);
    yuk_e_to_e(S_my, E[n1256], 0, -1, 2, scale);
    yuk_e_to_e(S_my, E[n12], 0, -1, 2, scale);
    yuk_e_to_e(S_py, E[sall], 0, 1, 2, scale);
    
    yuk_e_to_e(S_mx, E[eall], 0, 0, 2, scale);
    yuk_e_to_e(S_px, E[wall], 0, 0, 3, scale);
    yuk_e_to_e(S_px, E[w2468], 0, 0, 2, scale);
    yuk_e_to_e(S_px, E[w24], 0, 0, 2, scale);
    yuk_e_to_e(S_px, E[w2], 0, 0, 2, scale);
    break;
  case 2:
    yuk_e_to_e(S_mz, E[uall], 0, -1, 3, scale);
    yuk_e_to_e(S_mz, E[u1234], 0, -1, 2, scale);
    yuk_e_to_e(S_pz, E[dall], 0, 1, 2, scale);
    
    yuk_e_to_e(S_my, E[nall], 0, 0, 2, scale);
    yuk_e_to_e(S_py, E[sall], 0, 0, 3, scale);
    yuk_e_to_e(S_py, E[s3478], 0, 0, 2, scale);
    yuk_e_to_e(S_py, E[s34], 0, 0, 2, scale);
    
    yuk_e_to_e(S_mx, E[eall], 0, -1, 3, scale);
    yuk_e_to_e(S_mx, E[e1357], 0, -1, 2, scale);
    yuk_e_to_e(S_mx, E[e13], 0, -1, 2, scale);
    yuk_e_to_e(S_mx, E[e3], 0, -1, 2, scale);
    yuk_e_to_e(S_px, E[wall], 0, 1, 2, scale);
    break;
  case 3:
    yuk_e_to_e(S_mz, E[uall], -1, -1, 3, scale);
    yuk_e_to_e(S_mz, E[u1234], -1, -1, 2, scale);
    yuk_e_to_e(S_pz, E[dall], 1, 1, 2, scale);

    yuk_e_to_e(S_my, E[nall], 0, -1, 2, scale);
    yuk_e_to_e(S_py, E[sall], 0, 1, 3, scale);
    yuk_e_to_e(S_py, E[s3478], 0, 1, 2, scale);
    yuk_e_to_e(S_py, E[s34], 0, 1, 2, scale);
    
    yuk_e_to_e(S_mx, E[eall], 0, -1, 2, scale);
    yuk_e_to_e(S_px, E[wall], 0, 1, 3, scale);
    yuk_e_to_e(S_px, E[w2468], 0, 1, 2, scale);
    yuk_e_to_e(S_px, E[w24], 0, 1, 2, scale);
    yuk_e_to_e(S_px, E[w4], 0, 1, 2, scale);
    break;
  case 4:
    yuk_e_to_e(S_mz, E[uall], 0, 0, 2, scale);
    yuk_e_to_e(S_pz, E[dall], 0, 0, 3, scale);
    yuk_e_to_e(S_pz, E[d5678], 0, 0, 2, scale);
    
    yuk_e_to_e(S_my, E[nall], -1, 0, 3, scale);
    yuk_e_to_e(S_my, E[n1256], -1, 0, 2, scale);
    yuk_e_to_e(S_my, E[n56], -1, 0, 2, scale);
    yuk_e_to_e(S_py, E[sall], 1, 0, 2, scale);
    
    yuk_e_to_e(S_mx, E[eall], 1, 0, 3, scale);
    yuk_e_to_e(S_mx, E[e1357], 1, 0, 2, scale);
    yuk_e_to_e(S_mx, E[e57], 1, 0, 2, scale);
    yuk_e_to_e(S_mx, E[e5], 1, 0, 2, scale);
    yuk_e_to_e(S_px, E[wall], -1, 0, 2, scale);
    break;
  case 5:
    yuk_e_to_e(S_mz, E[uall], -1, 0, 2, scale);
    yuk_e_to_e(S_pz, E[dall], 1, 0, 3, scale);
    yuk_e_to_e(S_pz, E[d5678], 1, 0, 2, scale);
    
    yuk_e_to_e(S_my, E[nall], -1, -1, 3, scale);
    yuk_e_to_e(S_my, E[n1256], -1, -1, 2, scale);
    yuk_e_to_e(S_my, E[n56], -1, -1, 2, scale);
    yuk_e_to_e(S_py, E[sall], 1, 1, 2, scale);
    
    yuk_e_to_e(S_mx, E[eall], 1, 0, 2, scale);
    yuk_e_to_e(S_px, E[wall], -1, 0, 3, scale);
    yuk_e_to_e(S_px, E[w2468], -1, 0, 2, scale);
    yuk_e_to_e(S_px, E[w68], -1, 0, 2, scale);
    yuk_e_to_e(S_px, E[w6], -1, 0, 2, scale);
    break;
  case 6:
    yuk_e_to_e(S_mz, E[uall], 0, -1, 2, scale);
    yuk_e_to_e(S_pz, E[dall], 0, 1, 3, scale);
    yuk_e_to_e(S_pz, E[d5678], 0, 1, 2, scale);
    
    yuk_e_to_e(S_my, E[nall], -1, 0, 2, scale);
    yuk_e_to_e(S_py, E[sall], 1, 0, 3, scale);
    yuk_e_to_e(S_py, E[s3478], 1, 0, 2, scale);
    yuk_e_to_e(S_py, E[s78], 1, 0, 2, scale);
    
    yuk_e_to_e(S_mx, E[eall], 1, -1, 3, scale);
    yuk_e_to_e(S_mx, E[e1357], 1, -1, 2, scale);
    yuk_e_to_e(S_mx, E[e57], 1, -1, 2, scale);
    yuk_e_to_e(S_mx, E[e7], 1, -1, 2, scale);
    yuk_e_to_e(S_px, E[wall], -1, 1, 2, scale);
    break;
  case 7:
    yuk_e_to_e(S_mz, E[uall], -1, -1, 2, scale);
    yuk_e_to_e(S_pz, E[dall], 1, 1, 3, scale);
    yuk_e_to_e(S_pz, E[d5678], 1, 1, 2, scale);
    
    yuk_e_to_e(S_my, E[nall], -1, -1, 2, scale);
    yuk_e_to_e(S_py, E[sall], 1, 1, 3, scale);
    yuk_e_to_e(S_py, E[s3478], 1, 1, 2, scale);
    yuk_e_to_e(S_py, E[s78], 1, 1, 2, scale);
    
    yuk_e_to_e(S_mx, E[eall], 1, -1, 2, scale);
    yuk_e_to_e(S_px, E[wall], -1, 1, 3, scale);
    yuk_e_to_e(S_px, E[w2468], -1, 1, 2, scale);
    yuk_e_to_e(S_px, E[w68], -1, 1, 2, scale);
    yuk_e_to_e(S_px, E[w8], -1, 1, 2, scale);
    break;
  }

  yuk_e_to_l(S_mz, 'z', false, scale, L);
  yuk_e_to_l(S_pz, 'z', true, scale, L);
  yuk_e_to_l(S_my, 'y', false, scale, L);
  yuk_e_to_l(S_py, 'y', true, scale, L);
  yuk_e_to_l(S_mx, 'x', false, scale, L);
  yuk_e_to_l(S_px, 'x', true, scale, L);

  delete [] S;
}

void yuk_e_to_e(dcomplex_t *M, const dcomplex_t *W, int x, int y, int z,
                double scale) {
  const dcomplex_t *xs = builtin_yukawa_table_->xs(scale);
  const dcomplex_t *ys = builtin_yukawa_table_->ys(scale);
  const double *zs = builtin_yukawa_table_->zs(scale);
  const int *m = builtin_yukawa_table_->m(scale);
  const int *sm = builtin_yukawa_table_->sm(scale);
  int s = builtin_yukawa_table_->s();

  int offset = 0;
  for (int k = 0; k < s; ++k) {
    double factor_z = zs[4 * k + z];
    for (int j = 0; j < m[k] / 2; ++j) {
      int idx = (sm[k] + j) * 7 + 3;
      dcomplex_t factor_x = xs[idx + x];
      dcomplex_t factor_y = ys[idx + y];
      M[offset] += W[offset] * factor_z * factor_y * factor_x;
      offset++;
    }
  }
}

void yuk_e_to_l(const dcomplex_t *E, char dir, bool sgn, double scale, 
                dcomplex_t *L) {
  // Note: this function is called on the parent node.
  int p = builtin_yukawa_table_->p();
  int s = builtin_yukawa_table_->s();
  const int *M = builtin_yukawa_table_->m(scale);
  const int *sm = builtin_yukawa_table_->sm(scale);
  const int *smf = builtin_yukawa_table_->smf(scale);
  const int *f = builtin_yukawa_table_->f();
  const dcomplex_t *ealphaj = builtin_yukawa_table_->ealphaj(scale);
  double *legendre = new double[(p + 1) * (p + 2) / 2];
  const double *x = builtin_yukawa_table_->x();
  const double *w = builtin_yukawa_table_->w();
  double ld = builtin_yukawa_table_->lambda() *
    builtin_yukawa_table_->size(scale);
  const double *sqf = builtin_yukawa_table_->sqf();
  
  dcomplex_t *contrib = nullptr;
  dcomplex_t *W1 = new dcomplex_t[(p + 1) * (p + 2) / 2];
  dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];
  
  for (int k = 0; k < s; ++k) {
    int mk2 = M[k] / 2;
    
    // Compute sum_{j=1}^m(k) W(k, j) e^{-i * m * alpha_j}
    dcomplex_t *z = new dcomplex_t[f[k] + 1];
    
    // m = 0
    z[0] = 0.0;
    for (int j = 1; j <= mk2; ++j) {
      int idx = sm[k] + j - 1;
      z[0] += (E[idx] + conj(E[idx]));
    }
    
    // m = 1, ..., f[k], where m is odd
    for (int m = 1; m <= f[k]; m += 2) {
      z[m] = 0.0;
      for (int j = 1; j <= mk2; ++j) {
        int widx = sm[k] + j - 1;
        int aidx = smf[k] + (j - 1) * f[k] + m - 1;
        z[m] += (E[widx] - conj(E[widx])) * conj(ealphaj[aidx]);
      }
    }
    
    // m = 2, ..., f[k], where m is even
    for (int m = 2; m <= f[k]; m += 2) {
      z[m] = 0.0;
      for (int j = 1; j <= mk2; ++j) {
        int widx = sm[k] + j - 1;
        int aidx = smf[k] + (j - 1) * f[k] + m - 1;
        z[m] += (E[widx] + conj(E[widx])) * conj(ealphaj[aidx]);
      }
    }
    
    legendre_Plm_gt1_scaled(p, 1 + x[k] / ld, scale, legendre);
    
    double factor = w[k] / M[k];
    for (int n = 0; n <= p; ++n) {
      int mmax = (n <= f[k] ? n : f[k]);
      for (int m = 0; m <= mmax; ++m) {
        W1[midx(n, m)] += z[m] * legendre[midx(n, m)] * factor;
      }
    }
  }
  
  // Scale the local expansion by
  // (2 * n + 1) * (n - m)! / (n + m)! * (-1)^n * i^m * pi / 2 / ld
  int offset = 0;
  double factor = M_PI_2 / ld;
  for (int n = 0; n <= p; ++n) {
    dcomplex_t power_I{1.0, 0.0};
    for (int m = 0; m <= n; ++m) {
      W1[offset++] *= sqf[midx(n, m)] * power_I * factor;
      power_I *= dcomplex_t{0.0, 1.0};
    }
    factor *= -1;
  }
  
  if (!sgn) {
    // If the exponential expansion is not along the positive axis
    // direction with respect to the source, flip the sign of the
    // converted L_n^m where n is odd
    int offset = 1;
    for (int n = 1; n <= p; n += 2) {
      for (int m = 0; m <= n; ++m) {
        W1[offset++] *= -1;
      }
      offset += (n + 2);
    }
  }

  if (dir == 'z') {
    contrib = W1;
  } else if (dir == 'y') {
    const double *d = builtin_yukawa_table_->dmat_plus(0.0);
    yuk_rotate_sph_y(W1, d, W2);
    yuk_rotate_sph_z(W2, M_PI / 2, W1);
    contrib = W1;
  } else if (dir == 'x') {
    const double *d = builtin_yukawa_table_->dmat_minus(0.0);
    yuk_rotate_sph_y(W1, d, W2);
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
  delete [] legendre;
}


} //namespace dashmm

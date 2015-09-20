#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

#include "FMM_expansion.h"


namespace dashmm {

std::unique_ptr<Expansion> 
FMM_Expansion::S_to_M(Point center,
                      std::vector<Source>::iterator first,
                      std::vector<Source>::iterator last) const {
  FMM_Expansion *temp = new FMM_Expansion{center, p_, sqfr_, sqbinom_,
                                          dmat_plus_, dmat_minus_};

  std::complex<double> *expansion = temp->exp_data();
  const double *sqfr = sqfr_->data(); 
  double *legendre = new double[(p_ + 1) * (p_ + 2) / 2]; 

  for (auto i = first; i != last; ++i) {
    Point dist = i->position() - center;
    double q = i->charge();
    double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
    double r = dist.norm();

    double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

    // compute exp(-i * phi) for the azimuthal angle phi
    std::complex<double> ephi =
      (proj / r <= 1e-14 ? std::complex<double> (1.0, 0.0) :
       std::complex<double> (dist.x() / proj, -dist.y() / proj));

    // compute powers of r
    std::vector<double> powers_r(p_ + 1);
    powers_r[0] = 1.0;
    for (int j = 1; j <= p_; j++) {
      powers_r[j] = powers_r[j - 1] * r;
    }

    // compute powers of exp(-i * phi)
    std::vector<std::complex<double>> powers_ephi(p_ + 1);
    powers_ephi[0] = std::complex<double> (1.0, 0.0);
    for (int j = 1; j <= p_; j++) {
      powers_ephi[j] = powers_ephi[j - 1] * ephi;
    }

    // compute multipole expansion M_n^m
    legendre_Plm(p_, ctheta, legendre); 
    for (int m = 0; m <= p_; m++) {
      for (int n = m; n <= p_; n++) {
        expansion[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
          legendre[midx(n, m)] * sqfr[midx(n, m)]; 
      }
    }
  }

  delete [] legendre;

  return std::unique_ptr<Expansion>{temp};
}


std::unique_ptr<Expansion> 
FMM_Expansion::S_to_L(Point center,
                      std::vector<Source>::iterator first,
                      std::vector<Source>::iterator last) const {
  FMM_Expansion *temp = new FMM_Expansion{center, p_, sqfr_, sqbinom_,
                                          dmat_plus_, dmat_minus_};

  std::complex<double> *expansion = temp->exp_data();
  const double *sqfr = sqfr_->data(); 
  double *legendre = new double[(p_ + 1) * (p_ + 2) / 2];

  for (auto i = first; i != last; ++i) {
    Point dist = i->position() - center;
    double q = i->charge();
    double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
    double r = dist.norm();

    // compute cosine of the polar angle theta
    double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

    // compute exp(-i * phi) for the azimuthal angle phi
    std::complex<double> ephi =
      (proj / r <= 1e-14 ? std::complex<double>(1.0, 0.0) :
       std::complex<double> (dist.x() / proj, -dist.y() / proj));

    // compute powers of 1 / r
    std::vector<double> powers_r(p_ + 1);
    powers_r[0] = 1.0 / r;
    for (int j = 1; j <= p_; j++) {
      powers_r[j] = powers_r[j - 1] / r;
    }

    // compute powers of exp(-i * phi)
    std::vector<std::complex<double>> powers_ephi(p_ + 1);
    powers_ephi[0] = std::complex<double> (1.0, 0.0);
    for (int j = 1; j <= p_; j++) {
      powers_ephi[j] = powers_ephi[j - 1] * ephi;
    }

    // compute local expansion L_n^m
    legendre_Plm(p_, ctheta, legendre); 
    for (int m = 0; m <= p_; m++) {
      for (int n = m; n <= p_; n++) {
        expansion[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
          legendre[midx(n, m)] * sqfr[midx(n, m)]; 
      }
    }
  }
  delete [] legendre;

  return std::unique_ptr<Expansion>{temp};
}


std::unique_ptr<Expansion> 
FMM_Expansion::M_to_M(int from_child, double s_size) const {
  // The function is called on th expansion of the child box and
  // s_size is the child box's size.
  double h = s_size / 2;
  double px = center_.x() + (from_child % 2 == 0 ? h : -h);
  double py = center_.y() + (from_child % 4 <= 1 ? h : -h);
  double pz = center_.z() + (from_child < 4 ? h : -h);
  // The center of the returned expansion does not matter in the
  // following accumulation operation
  FMM_Expansion *temp = new FMM_Expansion{Point{px, py, pz}, p_,
                                          sqfr_, sqbinom_,
                                          dmat_plus_, dmat_minus_};

  // shift distance along the z-axis, combined with Y_n^0(pi, 0)
  const double rho = -sqrt(3) / 2 * s_size;

  // compute powers of rho
  std::vector<double> powers_rho(p_ + 1);
  powers_rho[0] = 1.0;
  for (int i = 1; i <= p_; i++) {
    powers_rho[i] = powers_rho[i - 1] * rho;
  }

  // allocate temp space to hold rotated spherical harmonic
  std::complex<double> *W1 = temp->exp_data();
  std::complex<double> *W2 = new std::complex<double>[(p_ + 2) * (p_ + 1) / 2];

  const std::complex<double> *MC = exp_.data();

  // table of rotation angle about the z-axis,
  // as an integer multiple of pi / 4
  const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
  // get rotation angle about the z-axis
  double alpha = tab_alpha[from_child] * M_PI_4;

  // get precomputed Wigner d-matrix for rotation about the y-axis
  const double *d1 = (from_child < 4 ?
                      dmat_plus_->at(1.0 / sqrt(3.0)).data() :
                      dmat_plus_->at(-1.0 / sqrt(3.0)).data());
  const double *d2 = (from_child < 4 ?
                      dmat_minus_->at(1.0 / sqrt(3.0)).data() :
                      dmat_minus_->at(-1.0 / sqrt(3.0)).data());

  const double *sqbinom = sqbinom_->data();

  // rotate the multipole expansion of the child box about z-axis
  rotate_sph_z(MC, alpha, p_, W1);

  // rotate the previous result further about the y-axis
  rotate_sph_y(W1, d1, p_, W2);

  // offset to operate multipole expansion
  int offset = 0;

  // shift along the z-axis by a distance of rho, write result in W1
  for (int n = 0; n <= p_; n++) {
    for (int m = 0; m <= n; m++) {
      W1[offset] = W2[offset];
      for (int k = 1; k <= n - m; k++) {
        W1[offset] += W2[midx(n - k, m)] * powers_rho[k] *
        sqbinom[midx(n - m, k)] * sqbinom[midx(n + m, k)];
      }
      offset++;
    }
  }

  // reverse rotate the shifted harmonic expansion about the y-axis
  rotate_sph_y(W1, d2, p_, W2);

  // reverse rotate the previous result further about the z-axis
  rotate_sph_z(W2, -alpha, p_, W1);

  delete [] W2;

  return std::unique_ptr<Expansion>{temp};
}


std::unique_ptr<Expansion> 
FMM_Expansion::M_to_L(Index s_index, double s_size, Index t_index) const {
  int t2s_x = s_index.x() - t_index.x();
  int t2s_y = s_index.y() - t_index.y();
  int t2s_z = s_index.z() - t_index.z();
  double tx = center_.x() - t2s_x * s_size;
  double ty = center_.y() - t2s_y * s_size;
  double tz = center_.z() - t2s_z * s_size;

  FMM_Expansion *temp = new FMM_Expansion{Point{tx, ty, tz}, p_,
                                          sqfr_, sqbinom_,
                                          dmat_plus_, dmat_minus_};
  // shifting distance
  double shift_unit = sqrt(t2s_x * t2s_x + t2s_y * t2s_y + t2s_z * t2s_z);
  const double rho = shift_unit * s_size;

  // compute powers of rho
  std::vector<double> powers_rho(p_ * 2 + 1);
  powers_rho[0] = 1.0 / rho;
  for (int i = 1; i <= p_ * 2; i++) {
    powers_rho[i] = powers_rho[i - 1] / rho;
  }


  // temp space to hold rotated spherical harmonic
  std::complex<double> *W1 = temp->exp_data();
  std::complex<double> *W2 = new std::complex<double>[(p_ + 2) * (p_ + 1) / 2];

  // get address of the source box's multipole expansion
  const std::complex<double> *M = exp_.data();

  // compute the projection of t2s on the x-y plane
  const double proj = sqrt(t2s_x * t2s_x + t2s_y * t2s_y);

  if (proj < 1e-14) {
    if (t2s_z > 0) {
      M_to_L_zp(M, powers_rho.data(), sqbinom_->data(), p_, W1);
    } else {
      M_to_L_zm(M, powers_rho.data(), sqbinom_->data(), p_, W1);
    }
  } else {
    // azimuthal angle
    double beta = acos(t2s_x / proj);
    if (t2s_y < 0) {
      beta = 2 * M_PI - beta;
    }

    // get precomputed Wigner d-matrix for rotation about y-axis
    const double *d1 = dmat_plus_->at(t2s_z / shift_unit).data();
    const double *d2 = dmat_minus_->at(t2s_z / shift_unit).data();

    rotate_sph_z(M, beta, p_, W1);
    rotate_sph_y(W1, d1, p_, W2);
    M_to_L_zp(W2, powers_rho.data(), sqbinom_->data(), p_, W1);
    rotate_sph_y(W1, d2, p_, W2);
    rotate_sph_z(W2, -beta, p_, W1);
  }

  delete [] W2;
  return std::unique_ptr<Expansion>{temp};
}


std::unique_ptr<Expansion> 
FMM_Expansion::L_to_L(int to_child, double t_size) const {
  // The function is called on the parent box and t_size is its child size
  double h = t_size / 2;
  double cx = center_.x() + (to_child % 2 == 0 ? -h : h);
  double cy = center_.y() + (to_child % 4 <= 1 ? -h : h);
  double cz = center_.z() + (to_child < 4 ? -h : h);

  // Again, the center of the expansion does not matter in the
  // following accumulation operation
  FMM_Expansion *temp = new FMM_Expansion{Point{cx, cy, cz}, p_,
                                          sqfr_, sqbinom_,
                                          dmat_plus_, dmat_minus_};
  // shift distance along the z-axis, combined with Y_n^0(pi, 0)
  const double rho = -sqrt(3) / 2 * t_size;

  // compute powers of rho
  std::vector<double> powers_rho(p_ + 1);
  powers_rho[0] = 1.0;
  for (int i = 1; i <= p_; i++) {
    powers_rho[i] = powers_rho[i - 1] * rho;
  }

  // temp space to hold rotated spherical harmonic
  std::complex<double> *W1 = temp->exp_data();
  std::complex<double> *W2 = new std::complex<double>[(p_ + 2) * (p_ + 1) / 2];

  // get address of the parent box's local expansion
  const std::complex<double> *LP = exp_.data();

  // table of rotation angle about the z-axis
  // as an integer multiple of pi / 4
  const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
  // get rotation angle about the z-axis
  double alpha = tab_alpha[to_child] * M_PI_4;

  // get precomputed Wigner d-matrix for rotation about the y-axis
  const double *d1 = (to_child < 4 ?
                      dmat_plus_->at(1.0 / sqrt(3.0)).data() :
                      dmat_plus_->at(-1.0 / sqrt(3.0)).data());
  const double *d2 = (to_child < 4 ?
                      dmat_minus_->at(1.0 / sqrt(3.0)).data() :
                      dmat_minus_->at(-1.0 / sqrt(3.0)).data());

  const double *sqbinom = sqbinom_->data();

  // offset to operate local expansion
  int offset = 0;

  // rotate the local expansion of the parent box about z-axis
  rotate_sph_z(LP, alpha, p_, W1);

  // rotate the previous result further about the y-axis
  rotate_sph_y(W1, d1, p_, W2);

  // shift along the z-axis by a distance of rho, write result in W1
  for (int n = 0; n <= p_; n++) {
    for (int m = 0; m <= n; m++) {
      W1[offset] = W2[offset];
      for (int k = 1; k <= p_ - n; k++) {
        W1[offset] += W2[midx(n + k, m)] * powers_rho[k] *
          sqbinom[midx(n + k - m, k)] * sqbinom[midx(n + k + m, k)];
      }
      offset++;
    }
  }

  // reverse rotate the shifted harmonic expansion about the y-axis
  rotate_sph_y(W1, d2, p_, W2);

  // reverse rotate the previous result further about the z-axis
  rotate_sph_z(W2, -alpha, p_, W1);

  delete [] W2;

  return std::unique_ptr<Expansion>{temp};
}


void FMM_Expansion::M_to_T(std::vector<Target>::iterator first,
                           std::vector<Target>::iterator last) const {
  const std::complex<double> *M = exp_.data();
  const double *sqfr = sqfr_->data(); 
  double *legendre = new double[(p_ + 1) * (p_ + 2) / 2];

  for (auto i = first; i != last; ++i) {
    Point dist = i->position() - center_;
    std::complex<double> potential = 0.0;
    double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
    double r = dist.norm();

    // compute cosine of the polar angle theta
    double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

    // compute exp(i * phi) for the azimuthal angle phi
    std::complex<double> ephi =
      (proj / r <= 1e-14 ? std::complex<double> (1.0, 0.0) :
       std::complex<double> (dist.x() / proj, dist.y() / proj));

    // compute powers of 1 / r
    std::vector<double> powers_r(p_ + 1);
    powers_r[0] = 1.0 / r;
    for (int j = 1; j <= p_; j++) {
      powers_r[j] = powers_r[j - 1] / r;
    }

    // compute powers of exp(i * phi)
    std::vector<std::complex<double>> powers_ephi(p_ + 1);
    powers_ephi[0] = std::complex<double>(1.0, 0.0);
    for (int j = 1; j <= p_; j++) {
      powers_ephi[j] = powers_ephi[j - 1] * ephi;
    }

    // evaluate the multipole expansion M_n^0
    legendre_Plm(p_, ctheta, legendre); 
    for (int n = 0; n <= p_; n++) {
      potential += M[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)]; 
    }

    // evaluate the multipole expansions M_n^m, where m = 1, ..., p
    for (int m = 1; m <= p_; m++) {
      for (int n = m; n <= p_; n++) {
        potential += 2.0 * real(M[midx(n, m)] * powers_ephi[m]) *
          powers_r[n] * legendre[midx(n, m)] * sqfr[midx(n, m)]; 
          //sqf[n - m] / sqf[n + m]; 
      }
    }

    i->set_phi(i->phi() + potential);
  }

  delete [] legendre;
}


void FMM_Expansion::L_to_T(std::vector<Target>::iterator first,
                           std::vector<Target>::iterator last) const {
  const std::complex<double> *L = exp_.data();
  const double *sqfr = sqfr_->data(); 
  double *legendre = new double[(p_ + 1) * (p_ + 2) / 2];

  for (auto i = first; i != last; ++i) {
    Point dist = i->position() - center_;
    std::complex<double> potential = 0.0;
    double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
    double r = dist.norm();

    // compute cosine of the polar angle theta
    double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

    // compute exp(i * phi) for the azimuthal angle phi
    std::complex<double> ephi =
      (proj / r <= 1e-14 ? std::complex<double>(1.0, 0.0) :
       std::complex<double>(dist.x() / proj, dist.y() / proj));

    // compute powers of r
    std::vector<double> powers_r(p_ + 1);
    powers_r[0] = 1.0;
    for (int j = 1; j <= p_; j++) {
      powers_r[j] = powers_r[j - 1] * r;
    }

    // compute powers of exp(i * phi)
    std::vector<std::complex<double>> powers_ephi(p_ + 1);
    powers_ephi[0] = std::complex<double>(1.0, 0.0);
    for (int j = 1; j <= p_; j++) {
      powers_ephi[j] = powers_ephi[j - 1] * ephi;
    }

    // evaluate the local expansion L_n^0
    legendre_Plm(p_, ctheta, legendre); 
    for (int n = 0; n <= p_; n++) {
      potential += L[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)]; 
    }

    // evaluate the local expansions L_n^m, where m = 1, ..., p
    for (int m = 1; m <= p_; m++) {
      for (int n = m; n <= p_; n++) {
        potential += 2.0 * real(L[midx(n, m)] * powers_ephi[m]) *
          powers_r[n] * legendre[midx(n, m)] * sqfr[midx(n, m)]; 
      }
    }

    i->set_phi(i->phi() + potential);
  }

  delete [] legendre;
}


void FMM_Expansion::S_to_T(std::vector<Source>::iterator s_first,
                           std::vector<Source>::iterator s_last,
                           std::vector<Target>::iterator t_first,
                           std::vector<Target>::iterator t_last) const {
  for (auto i = t_first; i != t_last; ++i) {
    std::complex<double> potential = 0.0;
    for (auto j = s_first; j != s_last; ++j) {
      Point s2t = i->position() - j->position();
      double dist = s2t.norm();
      if (dist > 0)
        potential += j->charge() / dist;
    }
    i->set_phi(i->phi() + potential);
  }
}


void FMM_Expansion::add_expansion(const Expansion *temp1) {
  for (size_t i = 0; i < temp1->size(); ++i)
    exp_[i] += temp1->term(i);
}


void FMM_Expansion::from_sum(const std::vector<const Expansion *> &exps) {
  //TODO:
}


std::unique_ptr<Expansion> FMM_Expansion::get_new_expansion(
    Point center) const {
  FMM_Expansion *temp = new FMM_Expansion{center, p_, sqfr_, sqbinom_,
                                          dmat_plus_, dmat_minus_};
  return std::unique_ptr<Expansion>{temp};
}

void FMM_Expansion::rotate_sph_z(const std::complex<double> *M,
                                 const double alpha, int p,
                                 std::complex<double> *MR) {
  // compute exp(i * alpha)
  std::complex<double> ealpha{cos(alpha), sin(alpha)};

  // compute powers of exp(i * alpha)
  auto powers_ealpha = new std::complex<double>[p + 1];
  powers_ealpha[0] = std::complex<double>(1.0, 0.0);
  for (int j = 1; j <= p; j++) {
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


void FMM_Expansion::rotate_sph_y(const std::complex<double> *M,
                                 const double *d, int p,
                                 std::complex<double> *MR) {
  int offset = 0;
  for (int n = 0; n <= p; n++) {
    int power_mp = 1;
    for (int mp = 0; mp <= n; mp++) {
      // retrieve address of wigner d-matrix entry d_n^{mp, 0}
      const double *coeff = &d[didx(n, mp, 0)];
      // get address of original harmonic expansion M_n^0
      const std::complex<double> *Mn = &M[midx(n, 0)];
      // compute rotated spherical harmonic M_n^mp
      MR[offset] = Mn[0] * coeff[0];
      double power_m = -1;
      for (int m = 1; m <= n; m++) {
        MR[offset] += (Mn[m] * power_m * coeff[m] + conj(Mn[m]) * coeff[-m]);
        power_m = -power_m;
      }
      MR[offset++] *= power_mp;
      power_mp = -power_mp;
    }
  }
}


void FMM_Expansion::M_to_L_zp(const std::complex<double> *M,
                              const double *rho, const double *sqbinom,
                              int p, std::complex<double> *L) {
  int offset = 0;
  for (int j = 0; j <= p; j++) {
    for (int k = 0; k <= j; k++) {
      L[offset] = 0;
      for (int n = k; n <= p; n++) {
        L[offset] += M[midx(n, k)] * pow_m1(n + k) * rho[j + n] *
          sqbinom[midx(n + j, n - k)] * sqbinom[midx(n + j, n + k)];
      }
      offset++;
    }
  }
}


void FMM_Expansion::M_to_L_zm(const std::complex<double> *M,
                              const double *rho, const double *sqbinom,
                              int p, std::complex<double> *L) {
  int offset = 0;
  for (int j = 0; j <= p; j++) {
    for (int k = 0; k <= j; k++) {
      L[offset] = 0;
      for (int n = k; n <= p; n++) {
        L[offset] += M[midx(n, k)] * pow_m1(k + j) * rho[j + n] *
          sqbinom[midx(n + j, n - k)] * sqbinom[midx(n + j, n + k)];
      }
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

std::vector<double> *generate_sqfr(const int p) {
  auto retval = new std::vector<double>((p + 1) * (p + 2) / 2); 

  double *temp = new double[2 * p + 1]; 
  temp[0] = 1.0; 
  for (int i = 1; i <= p * 2; i++) {
    temp[i] = temp[i - 1] * sqrt(i); 
  }

  double *sqfr = retval->data(); 
  for (int n = 0; n <= p; n++) {
    for (int m = 0; m <= n; m++) {
      sqfr[midx(n, m)] = temp[n - m] / temp[n + m]; 
    }
  }

  delete [] temp; 
  return retval;
}


std::vector<double> *generate_sqbinom(const int N) {
  int total = (N + 1) * (N + 2) / 2;
  auto retval = new std::vector<double>(total);
  double *val = retval->data();

  // mannually set the first three binomial coefficients
  val[0] = 1; // binom(0, 0)
  val[1] = 1; // binom(1, 0)
  val[2] = 1; // binom(1, 1)

  // compute the rest binomial coefficients
  for (int n = 2; n <= N; n++) {
    // get address of binom(n, 0) in current
    // get address of binom(n - 1, 0) in previous
    double *current = &val[midx(n, 0)];
    double *previous = &val[midx(n - 1, 0)];

    current[0] = 1; // binom(n, 0);
    for (int m = 1; m < n; m++) {
      current[m] = previous[m] + previous[m - 1];
    }
    current[n] = 1; // binom(n, n)
  }

  // compute the square root of the binomial coefficients
  for (int i = 0; i < total; i++) {
    val[i] = sqrt(val[i]);
  }

  return retval;
}


void generate_dmatrix_of_beta(const double beta, double *DP, double *DM,
                              const int p) {
  double cbeta = cos(beta);
  double sbeta = sin(beta);
  double s2beta2 = (1 - cbeta) / 2; // sin^2(beta / 2)
  double c2beta2 = (1 + cbeta) / 2; // cos^2(beta / 2)

  // set d_0^{0, 0} to 1
  DP[0] = 1;
  DM[0] = 1;

  // set d_1^{0, m}
  DP[1] = -sbeta / sqrt(2); // d_1^{0, -1}
  DP[2] = cbeta; // d_1^{0, 0}
  DP[3] = sbeta / sqrt(2); // d_1^{0, 1}
  DM[1] = -DP[1];
  DM[2] = DP[2];
  DM[3] = -DP[3];

  // set d_1^{1, m}
  DP[4] = s2beta2; // d_1^{1, -1}
  DP[5] = -sbeta / sqrt(2); // d_1^{1, 0}
  DP[6] = c2beta2; // d_1^{1, 1}
  DM[4] = DP[4];
  DM[5] = -DP[5];
  DM[6] = DP[6];

  // compute d_n^{0, m} for 2 <= n <= P
  for (int n = 2; n <= p; n++) {
    double *DPC = NULL, *DPP = NULL, *DMC = NULL;
    int m;
    // get address of d_n^{0, 0} for angle beta, saved in DPC
    // get address of d_{n - 1}^{0, 0} for angle beta, saved in DPP
    // get address of d_n^{0, 0} for angle -beta, saved in DMC
    DPC = &DP[didx(n, 0, 0)];
    DPP = &DP[didx(n - 1, 0, 0)];
    DMC = &DM[didx(n, 0, 0)];

    // compute d_n^{0, -n}
    m = -n;
    DPC[m] = -sbeta * 0.5 * sqrt((n - m) * (n - m - 1)) / n * DPP[m + 1];
    DMC[m] = DPC[m] * pow_m1(m);

    // compute d_n^{0, -n + 1}
    m = -n + 1;
    DPC[m] = (-sbeta * 0.5 * sqrt((n - m) * (n - m - 1)) * DPP[m + 1] +
              cbeta * sqrt((n - m) * (n + m)) * DPP[m]) / n;
    DMC[m] = DPC[m] * pow_m1(m);

    // handle d_n^{0, m}
    for (m = -n + 2; m <= n - 2; m++) {
      DPC[m] = (-sbeta * 0.5 * sqrt((n - m) * (n - m - 1)) * DPP[m + 1] +
                cbeta * sqrt((n - m) * (n + m)) * DPP[m] +
                sbeta * 0.5 * sqrt((n + m) * (n + m - 1)) * DPP[m - 1]) / n;
      DMC[m] = DPC[m] * pow_m1(m);
    }

    // handle d_n^{0, n - 1}
    m = n - 1;
    DPC[m] = (cbeta * sqrt((n - m) * (n + m)) * DPP[m] +
              sbeta * 0.5 * sqrt((n + m) * (n + m - 1)) * DPP[m - 1]) / n;
    DMC[m] = DPC[m] * pow_m1(m);

    // handle d_n^{0, n}
    m = n;
    DPC[m] = sbeta * 0.5 * sqrt((n + m) * (n + m - 1)) / n * DPP[m - 1];
    DMC[m] = DPC[m] * pow_m1(m);

    // compute d_n^{mp, m} for 1 <= mp <= n
    for (int mp = 1; mp <= n; mp++) {
      // get address of d_n^{mp, 0} for angle beta, saved in DPC
      // get address of d_{n - 1}^{mp - 1, 0} for angle beta, saved in DPP
      // get address of d_n^{mp, 0} for angle -beta, saved in DMC
      DPC = &DP[didx(n, mp, 0)];
      DPP = &DP[didx(n - 1, mp - 1, 0)];
      DMC = &DM[didx(n, mp, 0)];

      double factor = 1.0 / sqrt((n + mp) * (n + mp - 1));

      // compute d_n^{mp, -n}
      m = -n;
      DPC[m] = s2beta2 * sqrt((n - m) * (n - m - 1)) * factor * DPP[m + 1];
      DMC[m] = DPC[m] * pow_m1(mp - m);

      // compute d_n^{mp, -n + 1}
      m = -n + 1;
      DPC[m] = (s2beta2 * sqrt((n - m) * (n - m - 1)) * DPP[m + 1] -
                sbeta * sqrt((n + m) * (n - m)) * DPP[m]) * factor;
      DMC[m] = DPC[m] * pow_m1(mp - m);

      // compute d_n^{mp, m}
      for (m = -n + 2; m <= n - 2; m++) {
        DPC[m] = (s2beta2 * sqrt((n - m) * (n - m - 1)) * DPP[m + 1]
                  + c2beta2 * sqrt((n + m) * (n + m - 1)) * DPP[m - 1]
                  - sbeta * sqrt((n + m) * (n - m)) * DPP[m]) * factor;
        DMC[m] = DPC[m] * pow_m1(mp - m);
      }

      // compute d_n^{mp, n - 1}
      m = n - 1;
      DPC[m] = (-sbeta * sqrt((n + m) * (n - m)) * DPP[m] +
                c2beta2 * sqrt((n + m) * (n + m - 1)) * DPP[m - 1]) * factor;
      DMC[m] = DPC[m] * pow_m1(mp - m);

      // compute d_n^{mp, n}
      m = n;
      DPC[m] = c2beta2 * sqrt((n + m) * (n + m - 1)) * factor * DPP[m - 1];
      DMC[m] = DPC[m] * pow_m1(mp - m);
    }
  }
}


void generate_wigner_dmatrix(
    std::map<double, std::vector<double>, fmmtbl_cmp> *DP,
    std::map<double, std::vector<double>, fmmtbl_cmp> *DM, const int p) {
  double cbeta[24] = {1.0 / sqrt(5.0), 1.0 / sqrt(6.0), 1.0 / sqrt(9.0),
                      1.0 / sqrt(10.0), 1.0 / sqrt(11.0), 1.0 / sqrt(14.0),
                      1.0 / sqrt(19.0), 2.0 / sqrt(5.0), sqrt(2.0 / 3.0),
                      sqrt(1.0 / 2.0), sqrt(4.0 / 9.0), sqrt(1.0 / 3.0),
                      sqrt(4.0 / 13.0), sqrt(2.0 / 7.0), sqrt(4.0 / 17.0),
                      sqrt(2.0 / 11.0), sqrt(9.0 / 10.0), sqrt(9.0 / 11.0),
                      sqrt(9.0 / 13.0), sqrt(9.0 / 14.0), sqrt(9.0 / 17.0),
                      sqrt(9.0 / 19.0), sqrt(9.0 / 22.0), 0.0};

  int nd = (p + 1) * (4 * p * p + 11 * p + 6) / 6;
  for (int i = 0; i < 24; i++) {
    double beta = acos(cbeta[i]);
    std::vector<double> dp(nd), dm(nd);
    generate_dmatrix_of_beta(beta, dp.data(), dm.data(), p);
    (*DP)[cbeta[i]] = dp;
    (*DM)[cbeta[i]] = dm;
  }

  for (int i = 0; i < 23; i++) {
    double beta = acos(-cbeta[i]);
    std::vector<double> dp(nd), dm(nd);
    generate_dmatrix_of_beta(beta, dp.data(), dm.data(), p);
    (*DP)[-cbeta[i]] = dp;
    (*DM)[-cbeta[i]] = dm;
  }
}


} //namespace dashmm

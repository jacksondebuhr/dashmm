#ifndef __FMM_EXPANSION_CPP__
#define __FMM_EXPANSION_CPP__


#include <cmath>
#include <complex>
#include <map>
#include <vector>

#include "expansion.h"
#include "index.h"
#include "particle.h"
#include "point.h"


namespace dashmm {


struct fmmtbl_cmp {
  bool operator() (const double &a, const double &b) const {
    // The fmm table use cosine value as its key. The smallest gap
    // between key values is 0.01. Here, the code only compares the
    // first 6 digits and that should be enough.
    double aa = floor(a * 1000000.0) / 100000.0;
    double bb = floor(b * 1000000.0) / 100000.0;
    if (aa == bb)
      return false;
    return aa < bb;
  }
};


class FMM_Expansion : public Expansion {
 public:
  FMM_Expansion(Point center, int p, const std::vector<double> *sqfr,
      const std::vector<double> *sqbinom,
      const std::map<double, std::vector<double>, fmmtbl_cmp> *dmat_plus,
      const std::map<double, std::vector<double>, fmmtbl_cmp> *dmat_minus) :
    center_{center}, p_{p}, sqfr_{sqfr}, sqbinom_{sqbinom},
    dmat_plus_{dmat_plus}, dmat_minus_{dmat_minus} {
      exp_.resize((p + 2) * (p + 1) / 2, 0.0);
    }

  size_t size() const override {return exp_.size();}
  Point center() const override {return center_;}
  std::complex<double> term(size_t i) const override {return exp_[i];}

  std::unique_ptr<Expansion> S_to_M(Point center,
                                    std::vector<Source>::iterator first,
                                    std::vector<Source>::iterator last)
    const override;

  std::unique_ptr<Expansion> S_to_L(Point center,
                                    std::vector<Source>::iterator first,
                                    std::vector<Source>::iterator last)
    const override;

  std::unique_ptr<Expansion> M_to_M(int from_child,
                                    double s_size) const override;
  std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                    Index t_index) const override;
  std::unique_ptr<Expansion> L_to_L(int to_child, double t_size) const override;

  void M_to_T(std::vector<Target>::iterator first,
              std::vector<Target>::iterator last) const override;
  void L_to_T(std::vector<Target>::iterator first,
              std::vector<Target>::iterator last) const override;
  void S_to_T(std::vector<Source>::iterator s_first,
              std::vector<Source>::iterator s_last,
              std::vector<Target>::iterator t_first,
              std::vector<Target>::iterator t_last) const override;

  void add_expansion(const Expansion *temp1) override;
  void from_sum(const std::vector<const Expansion *> &exps) override;

  std::unique_ptr<Expansion> get_new_expansion(Point center) const override;

  std::complex<double> *exp_data() {return exp_.data();}
  const std::vector<double> *sqfr() const {return sqfr_;}
  const std::vector<double> *sqbinom() const {return sqbinom_;}

  const std::map<double, std::vector<double>, fmmtbl_cmp> *
  dmat_plus() const {
    return dmat_plus_;
  }

  const std::map<double, std::vector<double>, fmmtbl_cmp> *
  dmat_minus() const {
    return dmat_minus_;
  }

 private:
  Point center_;
  int p_;
  std::vector<std::complex<double>> exp_; // expansion
  const std::vector<double> *sqfr_;
  const std::vector<double> *sqbinom_;
  const std::map<double, std::vector<double>, fmmtbl_cmp> *dmat_plus_;
  const std::map<double, std::vector<double>, fmmtbl_cmp> *dmat_minus_;

  //void report_gsl_error(const int status) const;
  static void rotate_sph_z(const std::complex<double> *M, double alpha,
                           int p, std::complex<double> *MR);

  static void rotate_sph_y(const std::complex<double> *M, const double *d,
                           int p, std::complex<double> *MR);

  static void M_to_L_zp(const std::complex<double> *M, const double *rho,
                        const double *sqbinom, int p,
                        std::complex<double> *L);

  static void M_to_L_zm(const std::complex<double> *M, const double *rho,
                        const double *sqbinom, int p,
                        std::complex<double> *L);
};


std::vector<double> *generate_sqfr(const int p);

std::vector<double> *generate_sqbinom(const int n);

void generate_dmat_of_beta(const double beta, double *DP, double *DM,
                           const int p);

void generate_wigner_dmatrix(
    std::map<double, std::vector<double>, fmmtbl_cmp> *DP,
    std::map<double, std::vector<double>, fmmtbl_cmp> *DM, const int p);

void legendre_Plm(int n, double x, double *); 

inline int midx(const int n, const int m) {
  return n * (n + 1) / 2 + m;
}

inline int didx(const int n, const int mp, const int m) {
  return n * (n + 1) * (4 * n - 1) / 6 + mp * (2 * n + 1) + n + m;
}

inline double pow_m1(const int m) {
  return (m % 2 ? -1.0 : 1.0);
}


} //namespace dashmm

#endif

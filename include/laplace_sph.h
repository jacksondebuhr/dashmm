#ifndef __DASHMM_LAPLACE_SPH_EXPANSION_H__
#define __DASHMM_LAPLACE_SPH_EXPANSION_H__

/// \file include/laplace_sph.h
/// \brief Declaration of LaplaceSPH

#include <cmath>
#include <complex>
#include <map>
#include <vector>

#include "include/expansion.h"
#include "include/ids.h"
#include "include/index.h"
#include "include/particle.h"
#include "include/point.h"

namespace dashmm {

struct laplace_cmp {
  bool operator()(const double &a, const double &b) const {
    // The smallest gap between key values in the Laplace rotation matrix map is
    // 0.01. The operator compares the first 6 digits and that should be
    // enough. 
    double aa = floor(a * 1000000.0) / 100000.0;
    double bb = floor(b * 1000000.0) / 100000.0;
    if (aa == bb)
      return false;
    return aa < bb;
  }
};

using laplace_map_t = std::map<double, double *, laplace_cmp>; 

class LaplaceSPHTable {
public:
  LaplaceSPHTable(int n_digits); 
  ~LaplaceSPHTable(); 
  
  int p() const {return p_;}
  const double *sqf() const {return sqf_;}
  const double *sqbinom() const {return sqbinom_;}
  const double *dmat_plus(double v) const {return dmat_plus_->at(v);}
  const double *dmat_minus(double v) const {return dmat_minus_->at(v);}
  
private:
  int p_; 
  double *sqf_;
  double *sqbinom_;
  laplace_map_t *dmat_plus_; 
  laplace_map_t *dmat_minus_; 
  
  double *generate_sqf(); 
  double *generate_sqbinom(); 
  void generate_wigner_dmatrix(laplace_map_t *&dp, laplace_map_t *&dm); 
  void generate_dmatrix_of_beta(double beta, double *dp, double *dm); 
};

using uLaplaceSPHTable = std::unique_ptr<LaplaceSPHTable>; 
extern std::map<int, uLaplaceSPHTable> builtin_laplace_table_; 

struct LaplaceSPHData {
  int reserved;
  int type; 
  int n_digits; 
  Point center; 
  dcomplex_t expansion[];
}; 

class LaplaceSPH : public Expansion {
public:
  LaplaceSPH(Point center, int n_digits);
  LaplaceSPH(LaplaceSPHData *ptr, size_t bytes, int n_digits) {}
  ~LaplaceSPH(); 

  void *release() override {
    LaplaceSPHData *retval = data_; 
    data_ = nullptr; 
    return retval; 
  }
                      
  size_t bytes() const override {return bytes_;}
  bool valid() const override {return data_ != nullptr;}
  int type() const override {return kExpansionLaplaceSPH;}
  int accuracy() const override {return data_->n_digits;}
  bool provides_L() const override {return true;}
  bool provides_exp() const override {return false;}
  size_t size() const override {
    uLaplaceSPHTable &table = builtin_laplace_table_.at(data_->n_digits); 
    int p = table->p(); 
    return (p + 1) * (p + 2) / 2; 
  }

  Point center() const override {
    assert(valid()); 
    return data_->center; 
  }
  dcomplex_t term(size_t i) const override {
    return data_->expansion[i]; 
  }

  std::unique_ptr<Expansion> S_to_M(Point center, Source *first, 
                                    Source *last, double scale) const override; 
  std::unique_ptr<Expansion> S_to_L(Point center, Source *first, 
                                    Source *last, double scale) const override;
  std::unique_ptr<Expansion> M_to_M(int from_child, 
                                    double s_size) const override; 
  std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size, 
                                    Index t_index) const override; 
  std::unique_ptr<Expansion> L_to_L(int to_child, double t_size) const override;

  void M_to_T(Target *first, Target *last, double scale) const override; 
  void L_to_T(Target *first, Target *last, double scale) const override;
  void S_to_T(Source *s_first, Source *s_last, 
              Target *t_first, Target *t_last) const override; 

  void add_expansion(const Expansion *temp1) override; 
  std::unique_ptr<Expansion> get_new_expansion(Point center) const override; 
  
private: 
  LaplaceSPHData *data_; 
  size_t bytes_; 

  void rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR); 
  void rotate_sph_y(const dcomplex_t *M, const double *d, dcomplex_t *MR); 
  void M_to_L_zp(const dcomplex_t *M, const double *rho, double scale, 
                 dcomplex_t *L);  
  void M_to_L_zm(const dcomplex_t *M, const double *rho, double scale, 
                 dcomplex_t *L); 

                      
};

void legendre_Plm(int n, double x, double *P); 

inline int midx(const int n, const int m) {
  return n * (n + 1) / 2 + m;
}

inline int didx(const int n, const int mp, const int m) {
  return n * (n + 1) * (4 * n - 1) / 6 + mp * (2 * n + 1) + n + m;
}

inline double pow_m1(const int m) {
  return (m % 2 ? -1.0 : 1.0);
}

} // namespace dashmm

#endif

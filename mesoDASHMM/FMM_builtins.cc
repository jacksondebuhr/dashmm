#include "FMM_builtins.h"

#include <map>
#include <memory>
#include <vector>

#include "expansion.h"
#include "FMM_expansion.h"
#include "FMM_method.h"
#include "method.h"


namespace dashmm {


//TODO: make these heap allocations
FMM_Method builtin_fmm_method_{};
std::map<int, std::unique_ptr<Expansion>> builtin_fmm_expansion_{};
using fmm_table_t =
  std::tuple<std::vector<double> *,
             std::vector<double> *,
             std::map<double, std::vector<double>, fmmtbl_cmp> *,
             std::map<double, std::vector<double>, fmmtbl_cmp> *>;
std::map<int, fmm_table_t> builtin_fmm_tables_{};


Method *fmm_method() {
  return static_cast<Method *>(&builtin_fmm_method_);
}


Expansion *fmm_expansion(const int p) {
  auto it1 = builtin_fmm_expansion_.find(p);
  if (it1 != builtin_fmm_expansion_.end()) {
    return it1->second.get();
  } else {
    std::vector<double> *sqfr;
    std::vector<double> *sqbinom;
    std::map<double, std::vector<double>, fmmtbl_cmp> *dmat_plus;
    std::map<double, std::vector<double>, fmmtbl_cmp> *dmat_minus;

    auto it2 = builtin_fmm_tables_.find(p);
    if (it2 != builtin_fmm_tables_.end()) {
      std::tie(sqfr, sqbinom, dmat_plus, dmat_minus) = it2->second;
    } else {
      sqfr = generate_sqfr(p);
      sqbinom = generate_sqbinom(p * 2);
      dmat_plus = new std::map<double, std::vector<double>, fmmtbl_cmp>;
      dmat_minus = new std::map<double, std::vector<double>, fmmtbl_cmp>;
      generate_wigner_dmatrix(dmat_plus, dmat_minus, p);
      fmm_table_t entry(std::make_tuple(sqfr, sqbinom, dmat_plus, dmat_minus));
      builtin_fmm_tables_[p] = entry;
    }

    Expansion *added = new FMM_Expansion{Point{0.0, 0.0, 0.0}, p,
                                         sqfr, sqbinom, dmat_plus, dmat_minus};

    builtin_fmm_expansion_[p] = std::unique_ptr<Expansion>{added};
    return added;
  }
}


void fmm_builtins_init() {
}


void fmm_builtins_fini() {
  for (auto it = builtin_fmm_tables_.begin();
       it != builtin_fmm_tables_.end(); ++it) {
    std::vector<double> *sqfr;
    std::vector<double> *sqbinom;
    std::map<double, std::vector<double>, fmmtbl_cmp> *dmat_plus;
    std::map<double, std::vector<double>, fmmtbl_cmp> *dmat_minus;
    std::tie(sqfr, sqbinom, dmat_plus, dmat_minus) = it->second;
    delete sqfr;
    delete sqbinom;
    delete dmat_plus;
    delete dmat_minus;
  }
}


} // namespace dashmm

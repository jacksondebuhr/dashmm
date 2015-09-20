#include "BH_builtins.h"

#include <cassert>
#include <memory>
#include <vector>

#include "BH_expansion.h"
#include "BH_method.h"
#include "expansion.h"
#include "method.h"


namespace dashmm {


std::vector<std::unique_ptr<BH_Method>> *builtin_bh_method_{nullptr};
BH_Expansion *builtin_bh_expansion_{nullptr};


// a factory for the BH method
Method *bh_method(double crit) {
  //NOTE: This could be potentially wasteful if the user requests a given
  // method multiple times. But until such time as we have to, this will be
  // the implementation.
  BH_Method *added = new BH_Method{crit};
  builtin_bh_method_->emplace_back(added);
  return static_cast<Method *>(added);
}


Expansion *bh_expansion() {
  return static_cast<Expansion *>(builtin_bh_expansion_);
}


void bh_builtins_init() {
  builtin_bh_method_ = new std::vector<std::unique_ptr<BH_Method>>{};
  assert(builtin_bh_method_);
  builtin_bh_expansion_ = new BH_Expansion{Point{0.0, 0.0, 0.0}};
  assert(builtin_bh_expansion_);
}


void bh_builtins_fini() {
  if (builtin_bh_method_) {
    delete builtin_bh_method_;
    builtin_bh_method_ = nullptr;
  }
  if (builtin_bh_expansion_) {
    delete builtin_bh_expansion_;
    builtin_bh_expansion_ = nullptr;
  }
}


} //namespace dashmm

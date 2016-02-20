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


/// \file include/laplace_sph.h
/// \brief Implementation of laplace_sph_precompute()


#include "include/laplace_sph.h"


namespace dashmm {


void laplace_sph_precompute(int n_digits) {
  // If we are not running on SMP, the following needs to be wrapped and
  // broadcasted inside an hpx_run
  LaplaceSPHTableIterator entry = builtin_laplace_table_.find(n_digits);
  if (entry == builtin_laplace_table_.end()) {
    builtin_laplace_table_[n_digits] =
      std::unique_ptr<LaplaceSPHTable>{new LaplaceSPHTable{n_digits}};
  }
}


} // namespace dashmm

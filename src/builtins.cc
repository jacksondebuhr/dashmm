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


/// \file src/builtins.cc
/// \brief Implementation of built-in sources and expansions


#include "include/builtins.h"

#include <hpx/hpx.h>

#include "include/bh_method.h"
#include "include/direct_method.h"
#include "include/fmm_method.h"
#include "include/expansion.h"
#include "include/laplace_com.h"
#include "include/laplace_sph.h"
#include "include/method.h"



namespace dashmm {


Method *bh_method(double theta) {
  Method *retval = new BH{theta};
  return retval;
}


Method *bh_method_create_handler(size_t size, MethodSerial *data) {
  assert(size == (sizeof(MethodSerial) + sizeof(double)));
  assert(data->size == sizeof(double));
  double *theta = static_cast<double *>(data->data);
  return bh_method(*theta);
}
HPX_ACTION(HPX_FUNCTION, 0,
           bh_method_create_action, bh_method_create_handler,
           HPX_SIZE_T, HPX_POINTER);


Method *direct_method() {
  Method *retval = new Direct{};
  return retval;
}


Method *direct_method_create_handler(size_t size, MethodSerial *data) {
  assert(size == sizeof(MethodSerial));
  return direct_method();
}
HPX_ACTION(HPX_FUNCTION, 0,
           direct_method_create_action, direct_method_create_handler,
           HPX_SIZE_T, HPX_POINTER);


Method *fmm_method() {
  Method *retval = new FMM{};
  return retval;
}


Method *fmm_method_create_handler(size_t size, MethodSerial *data) {
  assert(size == sizeof(MethodSerial));
  return fmm_method();
}
HPX_ACTION(HPX_FUNCTION, 0,
           fmm_method_create_action, fmm_method_create_handler,
           HPX_SIZE_T, HPX_POINTER);


void register_built_in_methods() {
  assert(kSuccess == register_method(kMethodBH, bh_method_create_action, 0));
  assert(kSuccess == register_method(
                        kMethodDirect, direct_method_create_action, 0));
  assert(kSuccess == register_method(kMethodFMM, fmm_method_create_action, 0));
}


Expansion *laplace_COM_expansion() {
  Expansion *retval = new LaplaceCOM{Point{0.0, 0.0, 0.0}};
  return retval;
}


Expansion *laplace_COM_create_handler(double x, double y, double z,
                                      int n_digits) {
  LaplaceCOM *retval = new LaplaceCOM{Point{x, y, z}};
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           laplace_COM_create_action, laplace_COM_create_handler,
           HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE, HPX_INT);


Expansion *laplace_COM_interpret_handler(void *data, size_t size,
                                         int n_digits) {
  LaplaceCOM *retval = new LaplaceCOM{(LaplaceCOMData *)data, size};
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           laplace_COM_interpret_action, laplace_COM_interpret_handler,
           HPX_POINTER, HPX_SIZE_T, HPX_INT);


Expansion *laplace_sph_expansion(int n_digits) {
  // If we are not running on SMP, the following needs to be wrapped and
  // broadcasted inside an hpx_run
  LaplaceSPHTableIterator entry = builtin_laplace_table_.find(n_digits);
  if (entry == builtin_laplace_table_.end()) {
    builtin_laplace_table_[n_digits] =
      std::unique_ptr<LaplaceSPHTable>{new LaplaceSPHTable{n_digits}};
  }

  Expansion *retval = new LaplaceSPH{Point{0.0, 0.0, 0.0}, n_digits};
  return retval;
}


Expansion *laplace_sph_create_handler(double x, double y, double z,
                                      int n_digits) {
  LaplaceSPH *retval = new LaplaceSPH{Point{x, y, z}, n_digits};
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           laplace_sph_create_action, laplace_sph_create_handler,
           HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE, HPX_INT);


Expansion *laplace_sph_interpret_handler(void *data, size_t size,
                                         int n_digits) {
  LaplaceSPH *retval = new LaplaceSPH{(LaplaceSPHData *)data, size, n_digits};
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           laplace_sph_interpret_action, laplace_sph_interpret_handler,
           HPX_POINTER, HPX_SIZE_T, HPX_INT);


void register_built_in_expansions() {
  assert(kSuccess == register_expansion(kExpansionLaplaceCOM,
                                        laplace_COM_create_action,
                                        laplace_COM_interpret_action, 0));
  assert(kSuccess == register_expansion(kExpansionLaplaceSPH,
                                        laplace_sph_create_action,
                                        laplace_sph_interpret_action, 0));
}


}  // namespace dashmm

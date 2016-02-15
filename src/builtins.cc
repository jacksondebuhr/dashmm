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


}  // namespace dashmm

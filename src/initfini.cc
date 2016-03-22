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


/// \file src/initfini.cc
/// \brief Implemention of DASHMM initialization and finalization


#include <hpx/hpx.h>

#include "dashmm/types.h"


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


ReturnCode init(int *argc, char ***argv) {
  if (HPX_SUCCESS != hpx_init(argc, argv)) {
    return kRuntimeError;
  }

  return kSuccess;
}


ReturnCode finalize() {
  //shutdown the runtime
  hpx_finalize();

  return kSuccess;
}


} // namespace dashmm

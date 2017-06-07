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
/// \brief Implemention of DASHMM initialization and finalization


#include <hpx/hpx.h>
#include <libhpx/libhpx.h>

#include "dashmm/types.h"


namespace dashmm {


ReturnCode init(int *argc, char ***argv) {
  if (!hpx_initialized()) {
    if (HPX_SUCCESS != hpx_init(argc, argv)) {
      return kRuntimeError;
    }
  }

#ifdef DASHMM_INSTRUMENTATION
  if (libhpx_inst_tracer_active()) {
    libhpx_inst_phase_end();
  }
#endif

  return kSuccess;
}


ReturnCode finalize(bool shutdown) {
  if (shutdown) {
    hpx_finalize();
  }

  return kSuccess;
}


} // namespace dashmm

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


#include "dashmm/spmdutils.h"


/// \file
/// \brief Implementation of broadcast utility routine


namespace dashmm {


/// Action implementing broadcast
int broadcast_value_handler(char *value, size_t size) {
  hpx_exit(size, value);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           broadcast_value_action, broadcast_value_handler,
           HPX_POINTER, HPX_SIZE_T);


} // dashmm

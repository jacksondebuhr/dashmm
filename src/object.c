// =============================================================================
//  DASHMM
//
//  Copyright (c) 2014 - 2015, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license.  See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
//
//  Authors:
//    Jackson DeBuhr, Indiana University <jdebuhr [at] indiana.edu>
// =============================================================================


#include "libdashmm/object.h"

//other C stuff

//other dependencies

#include "libdashmm/basic_types.h"


uint32_t dashmm_object_class(dashmm_object_t *object) {
  return (object->properties & DASHMM_CLASS_MASK);
}


bool dashmm_user_object(dashmm_object_t *object) {
  return (object->properties & DASHMM_PROPERTY_USER);
}


bool dashmm_system_object(dashmm_object_t *object) {
  return (object->properties & DASHMM_PROPERTY_SYSTEM);
}


bool dashmm_verify_user_object(dashmm_handle_t handle, uint32_t type) {
  dashmm_object_t local;
  hpx_gas_memget_sync(&local, handle, sizeof(local));
  bool retval = true;
  if (dashmm_object_class(&local) != type) {
    retval = false;
  }
  if (!dashmm_user_object(&local)) {
    retval = false;
  }
  return retval;
}

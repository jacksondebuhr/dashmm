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


#include "libdashmm/array.h"

//c library includes

#include "hpx/hpx.h"

#include "libdashmm/object.h"


// ==================================================================
// Utilities
// ==================================================================


hpx_addr_t dashmm_array_offset(dashmm_handle_t handle, 
    uint64_t record, size_t offset) {
  dashmm_array_t local;
  hpx_gas_memget_sync(&local, handle, sizeof(local));
  
  if (local.num_records < record) {
    return HPX_NULL;
  }
  
  int64_t bytes = record * local.record_size + offset;
  return hpx_addr_add(local.data, bytes, local.block_size);
}


// ==================================================================
// Advanced Interface to DASHMM
// ==================================================================


dashmm_handle_t dashmm_array_alloc_distrib(uint64_t records, 
    size_t record_size,
    dashmm_array_distrib_t) {
  //
}


// ==================================================================
// Basic Interface to DASHMM
// ==================================================================


dashmm_handle_t dashmm_array_alloc(uint64_t records, size_t record_size) {
  return dashmm_array_alloc_distrib(records, 
                                    record_size, DASHMM_DISTRIB_CYCLIC);
}


int dashmm_array_free(dashmm_handle_t handle) {
  if (!dashmm_verify_user_object(handle, DASHMM_CLASS_ARRAY)) {
    return DASHMM_DOMAIN_ERROR;
  }
  
  dashmm_array_t local;
  hpx_gas_memget_sync(&local, handle, sizeof(local));
  hpx_addr_t done = hpx_lco_and_new(2);
  if (done == HPX_NULL) {
    return DASHMM_RUNTIME_ERROR;
  }
  
  hpx_gas_free(local.data, done);
  hpx_gas_free(handle, done);
  hpx_lco_wait(done);
  hpx_lco_delete(done, HPX_NULL);
  
  return DASHMM_SUCCESS;
}


int dashmm_array_memput(dashmm_handle_t handle, 
    uint64_t record, 
    size_t offset,
    size_t length,
    void *data) {
  if (!dashmm_verify_user_object(handle, DASHMM_CLASS_ARRAY)) {
    return DASHMM_DOMAIN_ERROR;
  }
  
  hpx_addr_t target = dashmm_array_offset(handle, record, offset);
  if (target == HPX_NULL) {
    return DASHMM_DOMAIN_ERROR;
  }
  
  hpx_addr_t done = hpx_lco_and_new(1);
  if (done == HPX_NULL) {
    return DASHMM_RUNTIME_ERROR;
  }
  
  hpx_gas_memput(target, data, length, HPX_NULL, done);
  hpx_lco_wait(done);
  hpx_lco_delete(done, HPX_NULL);
  
  return DASHMM_SUCCESS;
}


int dashmm_array_memget(dashmm_handle_t handle,
    uint64_t record,
    size_t offset,
    size_t length,
    void *data) {
  if (!dashmm_verify_user_object(handle, DASHMM_CLASS_ARRAY)) {
    return DASHMM_DOMAIN_ERROR;
  }
  
  hpx_addr_t source = dashmm_array_offset(handle, record, offset);
  if (source == HPX_NULL) {
    return DASHMM_DOMAIN_ERROR;
  }
  
  hpx_gas_memget_sync(data, source, length);
  
  return DASHMM_SUCCESS;
}



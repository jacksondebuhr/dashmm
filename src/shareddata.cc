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


// C / C++ stuff

#include "dashmm/shareddata.h"

// other DASHMM stuff


namespace dashmm {


namespace {
  void perform_reset(hpx_addr_t base, void *data, size_t bytes) {
    int nranks = hpx_get_num_ranks();
    hpx_addr_t done = hpx_lco_and_new(nranks);
    assert(done != HPX_NULL);

    for (int i = 0; i < nranks; ++i) {
      hpx_addr_t offset = hpx_addr_add(base, bytes * i, bytes);
      hpx_gas_memput_lsync(offset, data, bytes, done);
    }

    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
  }
}


int shared_data_construct_handler(size_t bytes, hpx_addr_t *retval) {
  *retval = hpx_gas_alloc_cyclic(hpx_get_num_ranks(), bytes, 0);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           shared_data_construct_action, shared_data_construct_handler,
           HPX_SIZE_T, HPX_POINTER);


int shared_data_destroy_handler(hpx_addr_t data) {
  hpx_gas_free_sync(data);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           shared_data_destroy_action, shared_data_destroy_handler,
           HPX_ADDR);


int shared_data_external_reset_handler(void *data, size_t UNUSED) {
  perform_reset(data->base, data->data, data->bytes);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           shared_data_external_reset_action,
           shared_data_external_reset_handler,
           HPX_POINTER, HPX_SIZE_T);


int shared_data_internal_reset_handler(void *data, size_t UNUSED) {
  perform_reset(data->base, data->data, data->bytes);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           shared_data_internal_reset_action,
           shared_data_internal_reset_handler,
           HPX_POINTER, HPX_SIZE_T);


} // dashmm

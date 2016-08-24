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


/// \file src/array.cc
/// \brief Implementation of DASHMM Array actions.


#include <cassert>
#include <cstring>

#include <hpx/hpx.h>

#include "dashmm/array.h"
#include "dashmm/reductionops.h"


namespace dashmm {


// allocate the array meta data
int allocate_array_meta_handler(void *UNUSED, size_t UNWANTED) {
  ArrayMetaAllocRunReturn retval{HPX_NULL, HPX_NULL, kSuccess};

  int ranks = hpx_get_num_ranks();
  retval.meta = hpx_gas_alloc_cyclic(ranks, sizeof(ArrayMetaData), 0);
  if (retval.meta == HPX_NULL) {
    retval.code = kAllocationError;
    hpx_exit(sizeof(retval), &retval);
  }

  retval.reducer = hpx_lco_reduce_new(ranks, sizeof(size_t) * ranks,
                                          size_sum_ident, size_sum_op);
  if (retval.reducer == HPX_NULL) {
    hpx_gas_free_sync(retval.meta);
    retval.meta = HPX_NULL;
    retval.code = kAllocationError;
  }

  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           allocate_array_meta_action, allocate_array_meta_handler,
           HPX_POINTER, HPX_SIZE_T);


// the local part of the allocation work
int allocate_local_work_handler(hpx_addr_t data, hpx_addr_t reducer,
                                size_t record_size, size_t record_count) {
  int ranks = hpx_get_num_ranks();
  int my_rank = hpx_get_my_rank();
  size_t *contrib = new size_t [ranks];
  assert(contrib);
  for (int i = 0; i < ranks; ++i) {
    if (i >= my_rank) {
      contrib[i] = record_count;
    } else {
      contrib[i] = 0;
    }
  }

  hpx_lco_set(reducer, sizeof(size_t) * ranks, contrib, HPX_NULL, HPX_NULL);
  hpx_lco_get(reducer, sizeof(size_t) * ranks, contrib);

  // Pin local part of the meta data
  hpx_addr_t global = hpx_addr_add(data,
                                   hpx_get_my_rank() * sizeof(ArrayMetaData),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));
  local->local_count = record_count;
  local->total_count = contrib[ranks - 1];
  local->size = record_size;
  local->data = HPX_NULL;

  delete [] contrib;

  int retval{0};
  if (record_count) {
    local->data = hpx_gas_alloc_local(1, record_count * record_size, 0);
    if (local->data == HPX_NULL) {
      retval = kAllocationError;
    }
  }

  hpx_gas_unpin(global);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_local_work_action, allocate_local_work_handler,
           HPX_ADDR, HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T);


// delete the reducer
int allocate_array_destroy_reducer_handler(hpx_addr_t reducer) {
  hpx_lco_delete_sync(reducer);
  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_array_destroy_reducer_action,
           allocate_array_destroy_reducer_handler,
           HPX_ADDR);


// Delete the local portion of the array
int deallocate_array_local_handler(hpx_addr_t meta) {
  hpx_addr_t global = hpx_addr_add(meta,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  if (local->data != HPX_NULL) {
    hpx_gas_free_sync(local->data);
  }

  hpx_gas_unpin(global);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           deallocate_array_local_action, deallocate_array_local_handler,
           HPX_ADDR);


/// Action that deallocates an array object
///
/// This will delete the ArrayMetaData and the records themselves.
///
/// \param obj - the global address of the array's meta data
///
/// \returns - HPX_SUCCESS
int deallocate_array_handler(hpx_addr_t obj) {
  hpx_bcast_rsync(deallocate_array_local_action, &obj);
  hpx_gas_free_sync(obj);
  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, 0, deallocate_array_action, deallocate_array_handler,
           HPX_ADDR);


/// Put data into an array object
///
/// This will fill records in a global array with the provided data. This
/// action is not well behaved as regards bad input arguments. This will likely
/// be updated in the future.
///
/// \param obj - the global address of the ArrayMetaData
/// \param first - the first record to put data into
/// \param last - the last (exclusive) record to put data into
/// \param in_data - the data to copy into the records
///
/// \returns - HPX_SUCCESS
int array_put_handler(hpx_addr_t obj, size_t first, size_t last,
                      void *in_data) {
  int retval = kSuccess;

  hpx_addr_t global = hpx_addr_add(obj,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  if (!hpx_gas_try_pin(global, (void **)&local)) {
    retval = kRuntimeError;
    hpx_exit(sizeof(retval), &retval);
  }

  char *data{nullptr};
  if (local->data == HPX_NULL) {
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);
  }
  if (!hpx_gas_try_pin(local->data, (void **)&data)) {
    hpx_gas_unpin(global);
    retval = kRuntimeError;
    hpx_exit(sizeof(retval), &retval);
  }

  // Do some simple bounds checking
  if (last > local->local_count || last < first) {
    retval = kDomainError;
    hpx_gas_unpin(local->data);
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);
  }

  char *beginning = data + first * local->size;
  size_t copy_size = (last - first) * local->size;
  memcpy(beginning, in_data, copy_size);

  hpx_gas_unpin(local->data);
  hpx_gas_unpin(global);

  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, 0, array_put_action, array_put_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER);


/// Get data from an array object
///
/// This will read records from a global array into the provided buffer. This
/// action is not well behaved as regards bad input arguments. This will likely
/// be updated in the future.
///
/// \param obj - the global address of the ArrayMetaData
/// \param first - the first record to read
/// \param last - the last (exclusive) record to read
/// \param out_data - a buffer into which the read data is stored
///
/// \returns - HPX_SUCCESS
int array_get_handler(hpx_addr_t obj, size_t first, size_t last,
                      void *out_data) {
  int retval = kSuccess;

  hpx_addr_t global = hpx_addr_add(obj,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  if (!hpx_gas_try_pin(global, (void **)&local)) {
    retval = kRuntimeError;
    hpx_exit(sizeof(retval), &retval);;
  }

  char *data{nullptr};
  if (local->data == HPX_NULL) {
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);;
  }
  if (!hpx_gas_try_pin(local->data, (void **)&data)) {
    retval = kRuntimeError;
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);;
  }

  // Do some simple bounds checking
  if (last > local->local_count || last < first) {
    retval = kDomainError;
    hpx_gas_unpin(local->data);
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);;
  }

  char *beginning = data + first * local->size;
  size_t copy_size = (last - first) * local->size;
  memcpy(out_data, beginning, copy_size);

  hpx_gas_unpin(local->data);
  hpx_gas_unpin(global);

  hpx_exit(sizeof(retval), &retval);;
}
HPX_ACTION(HPX_DEFAULT, 0, array_get_action, array_get_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER);


int array_local_count_handler(hpx_addr_t data) {
  hpx_addr_t global = hpx_addr_add(data,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  size_t retval = local->local_count;

  hpx_gas_unpin(global);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_local_count_action, array_local_count_handler,
           HPX_ADDR);


int array_total_count_handler(hpx_addr_t data) {
  hpx_addr_t global = hpx_addr_add(data,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  size_t retval = local->total_count;

  hpx_gas_unpin(global);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_total_count_action, array_total_count_handler,
           HPX_ADDR);

} // namespace dashmm

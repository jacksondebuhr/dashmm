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


// TODO: This will likely be obviated by some developments in hpx itself.
// For now, this will be left as the solution.
//
// NOTE: This means only one Array operation at a time.
namespace {
  hpx_addr_t allocate_save_reducer = HPX_NULL;
  hpx_addr_t allocate_save_meta = HPX_NULL;

  // set these local values. This is a junk workaround
  int allocate_save_some_data_handler(hpx_addr_t meta, hpx_addr_t reducer) {
    allocate_save_meta = meta;
    allocate_save_reducer = reducer;
    return HPX_SUCCESS;
  }
  HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
             allocate_save_some_data, allocate_save_some_data_handler,
             HPX_ADDR, HPX_ADDR);
}


// allocate the array meta data
int allocate_array_meta_handler(int *err) {
  int ranks = hpx_get_num_ranks();

  hpx_addr_t meta = hpx_gas_alloc_cyclic(ranks, sizeof(ArrayMetaData), 0);
  if (meta == HPX_NULL) {
    *err = kAllocationError;
    hpx_exit(HPX_ERROR);
  }

  hpx_addr_t reducer = hpx_lco_reduce_new(ranks, sizeof(size_t) * ranks,
                                          size_sum_ident, size_sum_op);
  if (reducer == HPX_NULL) {
    hpx_gas_free_sync(meta);
    *err = kAllocationError;
    hpx_exit(HPX_ERROR);
  }

  hpx_bcast_rsync(allocate_save_some_data, &meta, &reducer);

  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_array_meta_action, allocate_array_meta_handler,
           HPX_POINTER);


// the local part of the allocation work
int allocate_local_work_handler(hpx_addr_t *data, size_t record_size,
                                size_t record_count, int *err) {
  *data = allocate_save_meta;

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

  hpx_lco_set(allocate_save_reducer, sizeof(size_t) * ranks, contrib,
              HPX_NULL, HPX_NULL);
  hpx_lco_get(allocate_save_reducer, sizeof(size_t) * ranks, contrib);

  // Pin local part of the meta data
  hpx_addr_t global = hpx_addr_add(allocate_save_meta,
                                   hpx_get_my_rank() * sizeof(ArrayMetaData),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));
  local->local_count = record_count;
  local->total_count = contrib[ranks - 1];
  if (my_rank) {
    local->offset = contrib[my_rank - 1];
  } else {
    local->offset = 0;
  }
  local->size = record_size;
  local->data = HPX_NULL;

  if (record_count) {
    local->data = hpx_gas_alloc_local(1, record_count * record_size, 0);
    if (local->data == HPX_NULL) {
      *err = kAllocationError;
      int value = HPX_ERROR;
      return HPX_THREAD_CONTINUE(value);
    }
  }

  hpx_gas_unpin(global);

  int value = HPX_SUCCESS;
  return HPX_THREAD_CONTINUE(value);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_local_work_action, allocate_local_work_handler,
           HPX_POINTER, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER);


// delete the reducer and clear out the references to said reducer
int allocate_array_destroy_reducer_handler(void *UNUSED, size_t UNWANTED) {
  // delete the reducer
  hpx_lco_delete_sync(allocate_save_reducer);

  // then go ahead and clear out the values saved at each locality for safety
  hpx_addr_t null = HPX_NULL;
  hpx_bcast_rsync(allocate_save_some_data, &null, &null);

  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           allocate_array_destroy_reducer_action,
           allocate_array_destroy_reducer_handler,
           HPX_POINTER, HPX_SIZE_T);


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
  hpx_exit(HPX_SUCCESS);
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
int array_put_handler(hpx_addr_t obj, size_t first, size_t last, int *err,
                      void *in_data) {
  *err = kSuccess;

  hpx_addr_t global = hpx_addr_add(obj,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  if (!hpx_gas_try_pin(global, (void **)&local)) {
    *err = kRuntimeError;
    int retval = HPX_ERROR;
    return HPX_THREAD_CONTINUE(retval);
  }

  char *data{nullptr};
  if (local->data == HPX_NULL) {
    hpx_gas_unpin(global);
    int retval = HPX_SUCCESS;
    return HPX_THREAD_CONTINUE(retval);
  }
  if (!hpx_gas_try_pin(local->data, (void **)&data)) {
    hpx_gas_unpin(global);
    *err = kRuntimeError;
    int retval = HPX_ERROR;
    return HPX_THREAD_CONTINUE(retval);
  }

  // Do some simple bounds checking
  if (last > local->local_count || last < first) {
    *err = kDomainError;
    hpx_gas_unpin(local->data);
    hpx_gas_unpin(global);
    int retval = HPX_ERROR;
    return HPX_THREAD_CONTINUE(retval);
  }

  char *beginning = data + first * local->size;
  size_t copy_size = (last - first) * local->size;
  memcpy(beginning, in_data, copy_size);

  hpx_gas_unpin(local->data);
  hpx_gas_unpin(global);

  int retval = HPX_SUCCESS;
  return HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, 0, array_put_action, array_put_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER, HPX_POINTER);


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
int array_get_handler(hpx_addr_t obj, size_t first, size_t last, int *err,
                      void *out_data) {
  *err = kSuccess;

  hpx_addr_t global = hpx_addr_add(obj,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  if (!hpx_gas_try_pin(global, (void **)&local)) {
    *err = kRuntimeError;
    int retval = HPX_ERROR;
    return HPX_THREAD_CONTINUE(retval);
  }

  char *data{nullptr};
  if (local->data == HPX_NULL) {
    hpx_gas_unpin(global);
    int retval = HPX_SUCCESS;
    return HPX_THREAD_CONTINUE(retval);
  }
  if (!hpx_gas_try_pin(local->data, (void **)&data)) {
    *err = kRuntimeError;
    hpx_gas_unpin(global);
    int retval = HPX_ERROR;
    return HPX_THREAD_CONTINUE(retval);
  }

  // Do some simple bounds checking
  if (last > local->local_count || last < first) {
    *err = kDomainError;
    hpx_gas_unpin(local->data);
    hpx_gas_unpin(global);
    int retval = HPX_ERROR;
    return HPX_THREAD_CONTINUE(retval);
  }

  char *beginning = data + first * local->size;
  size_t copy_size = (last - first) * local->size;
  memcpy(out_data, beginning, copy_size);

  hpx_gas_unpin(local->data);
  hpx_gas_unpin(global);

  int retval = HPX_SUCCESS;
  return HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, 0, array_get_action, array_get_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER, HPX_POINTER);


} // namespace dashmm

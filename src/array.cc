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


namespace dashmm {


/// Action to allocate an array.
///
/// This will allocate both the ArrayMetaData object and the array storage
/// itself in the global address space.
///
/// \param count - the number of records to allocate
/// \param size - the size in bytes of each record
/// \param obj [out] - the address of the global address of the meta data
/// \param err [out] - the address of an integer in which to report errors
///
/// \returns - HPX_SUCCESS
int allocate_array_handler(size_t count, size_t size, hpx_addr_t *obj,
                           int *err) {
  *err = kSuccess;

  // create the metadata object
  hpx_addr_t retval = hpx_gas_alloc_local_at_sync(1, sizeof(ArrayMetaData), 0,
                                                  HPX_HERE);
  *obj = retval;
  if (retval == HPX_NULL) {
    *err = kAllocationError;
    hpx_exit(HPX_ERROR);
  }

  // NOTE: this should always work as the array meta data is allocated here.
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(retval, (void **)&meta)) {
    *err = kRuntimeError;
    hpx_gas_free_sync(retval);
    *obj = HPX_NULL;
    hpx_exit(HPX_ERROR);
  }
  meta->count = count;
  meta->size = size;

  meta->data = hpx_gas_alloc_local_at_sync(1, count * size, 0, HPX_HERE);
  if (meta->data == HPX_NULL) {
    *err = kAllocationError;
    hpx_gas_unpin(retval);
    hpx_gas_free_sync(retval);
    *obj = HPX_NULL;
    hpx_exit(HPX_ERROR);
  }

  hpx_gas_unpin(retval);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0, allocate_array_action, allocate_array_handler,
           HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER, HPX_POINTER);


/// Action that deallocates an array object
///
/// This will delete the ArrayMetaData and the records themselves.
///
/// \param obj - the global address of the array's meta data
///
/// \returns - HPX_SUCCESS
int deallocate_array_handler(hpx_addr_t obj) {
  // NOTE: This should always work as this is called on the root locality,
  // which is where the meta data was allocated.
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(obj, (void **)&meta)) {
    hpx_exit(HPX_ERROR);
  }

  assert(meta->data != HPX_NULL);
  hpx_gas_free_sync(meta->data);
  hpx_gas_unpin(obj);

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

  // NOTE: This should always work as this is called at the root locality, which
  // is where the meta data was allocated.
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(obj, (void **)&meta)) {
    *err = kRuntimeError;
    hpx_exit(HPX_ERROR);
  }

  char *local;
  if (!hpx_gas_try_pin(meta->data, (void **)&local)) {
    hpx_gas_unpin(obj);
    *err = kRuntimeError;
    hpx_exit(HPX_ERROR);
  }

  // Do some simple bounds checking
  if (last > meta->count || last < first) {
    *err = kDomainError;
    hpx_gas_unpin(meta->data);
    hpx_gas_unpin(obj);
    hpx_exit(HPX_ERROR);
  }

  char *beginning = local + first * meta->size;
  size_t copy_size = (last - first) * meta->size;
  memcpy(beginning, in_data, copy_size);

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(obj);
  hpx_exit(HPX_SUCCESS);
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

  // NOTE: This will always work as the action is called on the root locality,
  // which is where the meta data was allocated.
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(obj, (void **)&meta)) {
    *err = kRuntimeError;
    hpx_exit(HPX_ERROR);
  }

  char *local{nullptr};
  if (!hpx_gas_try_pin(meta->data, (void **)&local)) {
    *err = kRuntimeError;
    hpx_gas_unpin(obj);
    hpx_exit(HPX_ERROR);
  }

  // Do some simple bounds checking
  if (last > meta->count || last < first) {
    *err = kDomainError;
    hpx_gas_unpin(meta->data);
    hpx_gas_unpin(obj);
    hpx_exit(HPX_ERROR);
  }

  char *beginning = local + first * meta->size;
  size_t copy_size = (last - first) * meta->size;
  memcpy(out_data, beginning, copy_size);

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(obj);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0, array_get_action, array_get_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER, HPX_POINTER);


} // namespace dashmm

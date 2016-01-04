/// \file src/array.cc
/// \brief Implementation of DASHMM array object.


#include <cassert>
#include <cstring>

#include "hpx/hpx.h"

#include "include/array.h"
#include "include/types.h"


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


/// Action to allocate an array.
///
/// This will allocate both the ArrayMetaData object and the array storage
/// itself in the global address space.
///
/// \param count - the number of records to allocate
/// \param size - the size in bytes of each record
/// \param obj [out] - the address of the global address of the meta data
///
/// \returns - HPX_SUCCESS
int allocate_array_handler(size_t count, size_t size, hpx_addr_t *obj,
                           int *err) {
  *err = kSuccess;

  //create the metadata object
  hpx_addr_t retval = hpx_gas_alloc_local_at_sync(1, sizeof(ArrayMetaData), 0,
                                                  HPX_HERE);
  *obj = retval;
  if (retval == HPX_NULL) {
    *err = kAllocationError;
    hpx_exit(HPX_ERROR);
  }

  //pin it
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(retval, (void **)&meta)) {
    *err = kRuntimeError;
    hpx_exit(HPX_ERROR);
  }
  meta->count = count;
  meta->size = size;

  meta->data = hpx_gas_alloc_local_at_sync(1, count * size, 0, HPX_HERE);
  assert(meta->data != HPX_NULL);

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

  //Do some simple bounds checking
  if (last <= meta->count || last >= first) {
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
/// \param in_data - a buffer into which the read data is stored
///
/// \returns - HPX_SUCCESS
int array_get_handler(hpx_addr_t obj, size_t first, size_t last, int *err,
                      void *out_data) {
  *err = kSuccess;

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

  //Do some simple bounds checking
  if (last <= meta->count || last >= first) {
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


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


ReturnCode allocate_array(size_t records, size_t size, ObjectHandle *obj) {
  int runcode;
  int *arg = &runcode;
  hpx_run(&allocate_array_action, &records, &size, &obj, &arg);
  return runcode;
}


ReturnCode deallocate_array(ObjectHandle obj) {
  if (HPX_SUCCESS == hpx_run(&deallocate_array_action, &obj)) {
    return kSuccess;
  } else {
    return kRuntimeError;
  }
}


ReturnCode array_put(ObjectHandle obj, size_t first, size_t last,
                     void *in_data) {
  int runcode;
  int *arg = &runcode;
  hpx_run(&array_put_action, &obj, &first, &last, &arg, &in_data);
  return runcode;
}


ReturnCode array_get(ObjectHandle obj, size_t first, size_t last,
                     void *out_data) {
  int runcode;
  int *arg = &runcode;
  hpx_run(&array_get_action, &obj, &first, &last, &arg, &out_data);
  return runcode;
}


} // namespace dashmm

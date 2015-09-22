#include <cassert>
#include <cstring>

#include "hpx/hpx.h"

#include "include/types.h"


/////////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////////


struct ArrayMetaData {
  hpx_addr_t data;
  size_t count;
  size_t size;
};


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int allocate_array_handler(size_t count, size_t size, hpx_addr_t *obj) {
  //create the metadata object
  hpx_addr_t retval = hpx_gas_alloc(1, sizeof(ArrayMetaData), 0,
                                    HPX_DIST_TYPE_LOCAL);
  *obj = retval;
  if (retval == HPX_NULL) {
    hpx_exit(HPX_ERROR);
  }

  //pin it
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(retval, (void **)&meta)) {
    hpx_exit(HPX_ERROR);
  }
  meta->count = count;
  meta->size = size;

  meta->data = hpx_gas_alloc(1, count * size, 0, HPX_DIST_TYPE_LOCAL);
  hpx_gas_unpin(retval);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0, allocate_array_action, allocate_array_handler,
           HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER);


int deallocate_array_handler(hpx_addr_t obj) {
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(obj, (void **)&meta)) {
    hpx_exit(HPX_ERROR);
  }

  hpx_gas_free_sync(meta->data);
  hpx_gas_unpin(obj);

  hpx_gas_free_sync(obj);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0, deallocate_array_action, deallocate_array_handler,
           HPX_ADDR);


int array_put_handler(hpx_addr_t obj, size_t first, size_t last,
                      void *in_data) {
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(obj, (void **)&meta)) {
    hpx_exit(HPX_ERROR);
  }

  char *local;
  if (!hpx_gas_try_pin(meta->data, (void **)&local)) {
    hpx_gas_unpin(obj);
    hpx_exit(HPX_ERROR);
  }

  //Do some simple bounds checking
  assert(last <= meta->count);
  assert(last >= first);

  char *beginning = local + first * meta->size;
  size_t copy_size = (last - first) * meta->size;
  memcpy(beginning, in_data, copy_size);

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(obj);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0, array_put_action, array_put_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER);


int array_get_handler(hpx_addr_t obj, size_t first, size_t last,
                      void *out_data) {
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(obj, (void **)&meta)) {
    hpx_exit(HPX_ERROR);
  }

  char *local{nullptr};
  if (!hpx_gas_try_pin(meta->data, (void **)&local)) {
    hpx_gas_unpin(obj);
    hpx_exit(HPX_ERROR);
  }

  //Do some simple bounds checking
  assert(last <= meta->count);
  assert(last >= first);

  char *beginning = local + first * meta->size;
  size_t copy_size = (last - first) * meta->size;
  memcpy(out_data, beginning, copy_size);

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(obj);
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0, array_get_action, array_get_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER);


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


int allocate_array(size_t records, size_t size, ObjectHandle *obj) {
  int runcode = hpx_run(allocate_array_action, &records, &size, &obj);
  if (runcode != HPX_SUCCESS) {
    return kRuntimeError;
  } else if (obj == HPX_NULL) {
    return kAllocationError;
  }
  return kSuccess;
}


int deallocate_array(ObjectHandle obj) {
  if (HPX_SUCCESS == hpx_run(deallocate_array_action, &obj)) {
    return kSuccess;
  } else {
    return kRuntimeError;
  }
}


int array_put(ObjectHandle obj, size_t first, size_t last, void *in_data) {
  if (HPX_SUCCESS == hpx_run(array_put_action, &obj, &first, &last,
                            &in_data)) {
    return kSuccess;
  } else {
    return kRuntimeError;
  }
}


int array_get(ObjectHandle obj, size_t first, size_t last, void *out_data) {
  if (HPX_SUCCESS == hpx_run(array_get_action, &obj, &first, &last,
                             &out_data)) {
    return kSuccess;
  } else {
    return kRuntimeError;
  }
}

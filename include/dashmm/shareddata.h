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

#ifndef __DASHMM_SHARED_DATA_H__
#define __DASHMM_SHARED_DATA_H__


// TODO: proper documentation


#include <cassert>
#include <cstring>

#include <hpx/hpx.h>


namespace dashmm {

// Action for creating shared data objects
extern hpx_action_t shared_data_construct_action;

// Action for destroying shared data objects
extern hpx_action_t shared_data_destroy_action;

// Action for synchronous not inside hpx reset
extern hpx_action_t shared_data_external_reset_action;

// Action for asynchronous inside hpx reset
extern hpx_action_t shared_data_internal_reset_action;

// Marshalled type for reset actions
typedef struct reset_args_t {
  hpx_addr_t base;
  size_t bytes;
  char data[];
} reset_args_t;


template <typename T>
class LocalData {
 public:
  LocalData() : data_{HPX_NULL}, local_{nullptr} { }

  explicit LocalData(hpx_addr_t data) : data_{data}, local_{nullptr} {
    assert(data_ != HPX_NULL);
    assert(hpx_gas_try_pin(data_, (void **)local_));
  }

  LocalData(const LocalData<T> &other) {
    data_ = other.data_;
    assert(hpx_gas_try_pin(data_, (void **)local_));
  }

  LocalData(LocalData<T> &&other) {
    // Save the other here
    data_ = other.data_;
    local_ = other.local_;
    // Make sure other will not unpin
    other.local_ = nullptr;
    other.data_ = HPX_NULL;
  }

  ~LocalData() {
    if (data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
    }
  }

  LocalData<T> &operator=(const LocalData<T> &other) {
    // First remove original reference
    if (data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
    }

    // Then create updated reference
    data_ = other.data_;
    assert(hpx_gas_try_pin(data_, (void **)local_));

    return *this;
  }

  LocalData<T> &operator=(LocalDat<T> &&other) {
    if (data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
    }

    data_ = other.data_;
    local_ = other.local_;

    other.data_ = HPX_NULL;
    other.local_ = nullptr;

    return *this;
  }

  const T& operator*() const {return *local_;}
  const T* operator->() const {return local_;}

 private:
   hpx_addr_t data_;
   T *local_;
};


// This object is used to represent read-only (at least during DASHMM library
// calls) data. The data is replicated at each locality.
//
// The intention of such a thing is to provide a way to have quick local access
// to read only data that is used frequently. Otherwise, if this were a normal
// global address space object, we would be using the network a great deal.
//
// The implementation will be as a cyclic allocation, with one block of the
// size of T for each locality. Then the value() member will just go ahead and
// compute the correct offset from the saved base, and we can pin and unpin
// as needed.
template <typename T>
class SharedData {
 public:
  SharedData(const T *value) : data_{HPX_NULL} {
    if (hpx_is_active()) {
      data_ = hpx_gas_alloc_cyclic(hpx_get_num_ranks(), sizeof(T), 0);
    } else {
      hpx_addr_t *argaddr = &data_;
      size_t bytes = sizeof(T);
      hpx_run(&shared_data_construct_action, &bytes, &argaddr);
    }
    assert(data_ != HPX_NULL);
  }

  SharedData(hpx_addr_t data) : data_{data} { }

  // We do not do this with a destructor becuase these objects might
  // exist at the end of their containing scope.
  //
  // The user must assure that all LocalData objects are destroyed prior to
  // this object, otherwise there will be a use-after-free error.
  void destroy() {
    if (hpx_is_active()) {
      hpx_gas_free_sync(data_);
    } else {
      assert(!hpx_is_active());   // It is an error to use this from inside HPX
      hpx_run(&shared_data_destroy_action, &data_);
    }
    data_ = HPX_NULL;
  }

  // Set the value - this is synchronous
  void reset(const T *value) {
    if (hpx_is_active()) {
      hpx_addr_t done = reset_async(value);
      hpx_lco_wait(done);
      hpx_lco_delete_sync(done);
    } else {
      size_t argsize = sizeof(reset_args_t) + sizeof(T);
      reset_args_t *parms = new char[argsize];
      parms.base = data_;
      parms.bytes = sizeof(T);
      memcpy(parms.data, value, sizeof(T));
      hpx_run(&shared_data_external_reset_action, parms, argsize);
      delete [] parms;
    }
  }

  // This returns an LCO that indicates completion of the local caching.
  // The caller must free that LCO. This is callable from HPX threads, but not
  // outside them.
  hxp_addr_t reset_async(const T *value) {
    assert(hpx_is_active());

    size_t argsize = sizeof(reset_args_t) + sizeof(T);
    reset_args_t *parms = new char[argsize];
    parms.base = data_;
    parms.bytes = sizeof(T);
    memcpy(parms.data, value, sizeof(T));

    hpx_addr_t retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);

    hpx_call(HPX_HERE, shared_data_internal_reset_action, retval,
             parms, argsize);

    delete [] parms;

    return retval;
  }

  // Get the local object - this should only be called from inside HPX threads
  LocalData<T> value() const {
    assert(hpx_is_active());
    hpx_addr_t offset = hpx_addr_add(data_, sizeof(T) * hpx_get_my_rank(),
                                     sizeof(T));
    return LocalData<T>{offset};
  }

  // Get the base of the array
  hpx_addr_t data() const {return data_;}

 private:
  hpx_addr_t data_;
};


} // dashmm


#endif // __DASHMM_SHARED_DATA_H__

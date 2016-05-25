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


/// \file include/dashmm/shareddata.h
/// \brief Definition of SharedData type.


#include <cassert>
#include <cstring>

#include <hpx/hpx.h>


namespace dashmm {


/// Action for creating shared data objects
extern hpx_action_t shared_data_construct_action;

/// Action for destroying shared data objects
extern hpx_action_t shared_data_destroy_action;

/// Action for synchronous not inside hpx reset
extern hpx_action_t shared_data_external_reset_action;

/// Action for asynchronous inside hpx reset
extern hpx_action_t shared_data_internal_reset_action;


/// Marshalled type for reset actions
struct reset_args_t {
  hpx_addr_t base;
  size_t bytes;
  char data[];
};


/// Local reference to a SharedData object
///
/// LocalData is the version of a SharedData object that can be used
/// on a particular locality. See the description of SharedData for more.
/// This object predominantly is tasked with serving the pinned address that
/// is relevant for the locality on which this object was created. This acts
/// like a pointer type, and provides both * and -> operations to make use
/// of the shared data easier.
///
/// This object captures resources (it pins a global address). When this
/// object is destroyed, the pin is removed.
///
/// This object is usable only from inside an HPX thread.
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

  LocalData<T> &operator=(LocalData<T> &&other) {
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

  /// Return the local address translation of the represented object.
  ///
  /// This is provided to avoid having to pin/unpin a great number of times
  /// if a given LocalData object needs to be passed to a number of functions,
  /// as occurs somewhere in the library. In practice, this should be avoided
  /// for simple uses.
  T *value() const {return local_;}

 private:
   hpx_addr_t data_;
   T *local_;
};


/// Represents data shared among the localities during an evaluation
///
/// This object serves the following case: data that is read frequently at
/// each locality, but is only rarely modified. Implemented as a normal
/// global address space object, this use case would require a good deal of
/// network traffic for data that does not change. Once set up, this object
/// will obviate the need for that traffic, and will instead provide a means to
/// get access to a local copy of the data relatively cheaply.
///
/// The template parameter gives the type that is shared around the system.
/// The only requirement on that data is that is be trivially copyable.
///
/// Local access can be obtained through the value() method.
///
/// This object can be used from both inside HPX threads and from outside
/// HPX threads. With two exceptions, mentioned below, that is true of each
/// method assocaited with this class.
template <typename T>
class SharedData {
 public:
  /// Construct from a value
  ///
  /// This will allocate the global memory and associates it with this object.
  /// If the provided data is a null pointer, no value will be shared across
  /// the localities, though the global memory will have been allocated.
  ///
  /// \param value - the data to share across the system.
  SharedData(const T *value) : data_{HPX_NULL} {
    if (hpx_is_active()) {
      data_ = hpx_gas_alloc_cyclic(hpx_get_num_ranks(), sizeof(T), 0);
    } else {
      hpx_addr_t *argaddr = &data_;
      size_t bytes = sizeof(T);
      hpx_run(&shared_data_construct_action, &bytes, &argaddr);
    }
    assert(data_ != HPX_NULL);
    if (value != nullptr) {
      reset(value);
    }
  }

  /// Construct from an HPX address.
  ///
  /// This does not allocate new memory. Instead, this is used to interpret
  /// a given global address as a SharedData object. This cannot assure that
  /// the given address actually contains data for an object of type T.
  SharedData(hpx_addr_t data) : data_{data} { }

  /// Destroy the global data backing this object.
  ///
  /// The lifetime of SharedData objects is expected to be longer than their
  /// enclosing scope, so we do not destroy the global memory when this
  /// object is destroyed. Instead, the user must call destroy() explicitly.
  ///
  /// It is an error to have a LocalData object associated with object when
  /// destroy() is called. Having such an outstanding LocalData object will
  /// lead to use-after-free errors. It is the user's responsibility to assure
  /// that this does not happen.
  ///
  /// One means to assure that the above does not occur is to only call this
  /// method from outside an HPX-5 thread. LocalData objects can only be used
  /// from inside HPX threads, and so being outside any HPX thread will assure
  /// that all LocalData objects associated with this have been destroyed.
  ///
  /// That being said, it is possible to use this routine from inside HPX-5,
  /// provided the user remains aware of the previous warning.
  void destroy() {
    if (hpx_is_active()) {
      hpx_gas_free_sync(data_);
    } else {
      assert(!hpx_is_active());   // It is an error to use this from inside HPX
      hpx_run(&shared_data_destroy_action, &data_);
    }
    data_ = HPX_NULL;
  }

  /// Reset the value of the shared data
  ///
  /// This will reset the value shared around the system to that provded
  /// by the argument. This method is synchronous: it will not return until
  /// the data has been updated at each locality.
  ///
  /// \param value - the address of the object to share. Cannot be null.
  void reset(const T *value) {
    assert(value != nullptr);
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

  /// Reset the value of the shared data.
  ///
  /// This will reset the valuye shared around the system to that provided
  /// by the arguments. This method is asynchronous: it can return before the
  /// data has been updated everywhere. The returned address is for an LCO that
  /// can be used to wait on completion of the update.
  ///
  /// This cannot be called from outside an HPX-5 thread. The caller assumes
  /// ownership of the returned LCO, and the called is responsible for freeing
  /// the LCO when it is no longer needed.
  ///
  /// \param value - the addess of the object to share. Cannot be null.
  ///
  /// \returns - address of an LCO that indicates completion of the reset.
  hxp_addr_t reset_async(const T *value) {
    assert(hpx_is_active());
    assert(value != nullptr);

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

  /// Create a local access object for the calling locality
  ///
  /// This will return a LocalData object that refers to the shared data for
  /// the calling locality.
  ///
  /// This can only be used from inside an HPX-5 thread.
  LocalData<T> value() const {
    assert(hpx_is_active());
    hpx_addr_t offset = hpx_addr_add(data_, sizeof(T) * hpx_get_my_rank(),
                                     sizeof(T));
    return LocalData<T>{offset};
  }

  /// Return the global address of the backing memory
  hpx_addr_t data() const {return data_;}

 private:
  hpx_addr_t data_;
};


} // dashmm


#endif // __DASHMM_SHARED_DATA_H__

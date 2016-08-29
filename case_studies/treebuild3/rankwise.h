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

#ifndef __DASHMM_RANKWISE_H__
#define __DASHMM_RANKWISE_H__


#include <hpx/hpx.h>


/// This object represents the local portion of a RankWise object.
///
/// This object will represent either the data on the rank it is being
/// used, or it will pull a copy across the network if the object is created
/// for a different rank.
template <typename T>
class RankLocal {
public:
  /// Construct a RankLocal object from a given address.
  ///
  /// This will attempt to pin the data. If that fails, the assumption is made
  /// is that the data is remote, and so a memget is issued to retrieve the
  /// information.
  ///
  /// \param - global address of data of type T
  RankLocal(hpx_addr_t global)
        : local_{nullptr}, global_{global}, remote_{false} {
    if (global_ == HPX_NULL) return;

    if (!hpx_gas_try_pin(global_, (void **)&local_)) {
      //If the pin fails, assume that it is because the address is not local
      // so get that data
      local_ = static_cast<T *>(hpx_malloc_registered(sizeof(T)));
      if (local_) {
        remote_ = true;
        hpx_gas_memget_sync(local_, global_, sizeof(T));
      } else {
        global_ = HPX_NULL;
      }
    }
  }

  /// Release any acquired resources.
  ~RankLocal() {
    if (local_ != nullptr) {
      if (remote_) {
        hpx_free_registered(local_);
      } else {
        hpx_gas_unpin(global_);
      }
      global_ = HPX_NULL;
      local_ = nullptr;
    }
  }

  /// We delete copy assignment and construction to simplify reference counting.
  RankLocal(const RankLocal &other) = delete;
  RankLocal &operator=(const RankLocal &other) = delete;

  /// Moves are allowed. The target of the move assumes all resources.
  RankLocal(RankLocal &&other) {
    local_ = other.local_;
    global_ = other.global_;
    remote_ = other.remote_;
    other.local_ = nullptr;
    other.global_ = HPX_NULL;
    other.remote_ = false;
  }

  /// Moves are allowed. The target of the move assumes all resources.
  RankLocal &operator=(RankLocal &&other) {
    local_ = other.local_;
    global_ = other.global_;
    remote_ = other.remote_;
    other.local_ = nullptr;
    other.global_ = HPX_NULL;
    other.remote_ = false;
  }

  /// A RankLocal object is valid if there is a non null global address,
  /// and a local version was acquired in some fashion.
  bool valid() const {return global_ != HPX_NULL && local_ != nullptr;}
  bool remote() const {return remote_;}

  /// Typical dereferencing operations to allow for use of the RankLocal as
  /// if it were a pointer to the underlying type.
  T& operator*() {return *local_;}
  T* operator->() {return local_;}

private:
  T *local_;
  hpx_addr_t global_;
  bool remote_;
};


/// Rankwise data
///
/// This object represents the concept of needing a particular set of data
/// at each rank, but which might need to point to different addresses on each
/// rank. The prototypical example of this will be the tree data: each rank
/// builds a local tree, but each local tree is representing the same
/// hierarchical subdivision, the local trees are bound together into a single
/// RankWise object.
///
/// The implementation is via a cyclic allocation that can store a particular
/// type of data. With the global address of the RankWise allocation, any
/// rank can get to the data at any other locality.
template <typename T>
class RankWise {
 public:
  /// Construct a RankWise object.
  ///
  /// This will either create an empty RankWise object (which must then have
  /// allocate() called on it) or will associate this object with the given
  /// global address.
  ///
  /// \param dat - an address to associate with this object, can be HPX_NULL
  RankWise(hpx_addr_t dat = HPX_NULL) : data_{dat} { }

  /// Allocate a RankWise object
  ///
  /// A default constructed RankWise object is not associated with anything
  /// in the global address space. This will allocate memory with the correct
  /// size and placement to be used as a RankWise object.
  ///
  /// This initializes the memory to all zero values. The user is required to
  /// set other initial values if that is needed.
  void allocate() {
    assert(data_ == HPX_NULL);
    data_ = hpx_gas_calloc_cyclic(hpx_get_num_ranks(), sizeof(T), 0);
  }

  /// Is the RankWise object valid
  ///
  /// A RankWise object is valid if it is associated with some piece of global
  /// memory. Of course, this does not test that the memory referred actually
  /// serves as a RankWise object, but without much more involved testing that
  /// would potentially use the network, this is a good first approximation.
  ///
  /// \returns - true if the object is valid; false otherwise
  bool valid() const {return data_ != HPX_NULL;}

  /// Return the global address represented by this object
  hpx_addr_t data() const {return data_;}

  /// Return the local version of the represented data
  ///
  /// This returns a RankLocal object representing the calling ranks portion
  /// of the RankWise data. See RankLocal for details.
  RankLocal<T> here() const {
    return there(hpx_get_my_rank());
  }

  /// Return the represented data for a given rank.
  ///
  /// This returns a RankLocal object representing the given rank's portion
  /// of the RankWise data. See RankLocal for details.
  RankLocal<T> there(int rank) const {
    return RankLocal<T>(address_at(rank));
  }

 private:
  /// Utility member to perform address arithmetic.
  hpx_addr_t address_at(int rank) const {
    assert(data_ != HPX_NULL);
    return hpx_addr_add(data_, sizeof(T) * rank, sizeof(T));
  }

  /// The global address of the data serving the RankWise object
  hpx_addr_t data_;
};


#endif // __DASHMM_RANKWISE_H__

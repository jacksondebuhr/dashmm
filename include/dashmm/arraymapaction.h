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


#ifndef __DASHMM_ARRAY_MAP_ACTION_H__
#define __DASHMM_ARRAY_MAP_ACTION_H__


/// \file include/arraymapaction.h
/// \brief Definitions for actions mappable over an array.


#include <hpx/hpx.h>

#include "dashmm/arraymetadata.h"


namespace dashmm {


// Forward declare Array so we can be friends with it.
template <typename T>
class Array;


/// ArrayMapAction
///
/// A class that represents an action that can be invoked on the elements of
/// a DASHMM array. It is a template requiring both the type for the array
/// records as well as an environment that is passed to each invocation of the
/// associated action.
///
/// To specify the action, the user will provide a function pointer with the
/// following signature: void func(T *, const size_t, const size_t, const E *),
/// which has been aliased as map_function_t. When writing such a function,
/// it is important to know that the meaning of the provided arguments:
///
/// * The first is a pointer to an array of T objects. These are the data on
///   which the action will act.
/// * The second is the total count of records to be examined in the current
///   call of the function.
/// * The third is the offset of the provided pointer in the overall array.
///   NOTE: This does not mean that one should begin indexing the provided
///   array at the value of the third argument. Instead, this is provided in
///   case there is some reason that it is important to know the overall offset
///   in the array.
/// * The final is the environment passed into the Array<T>::map() function
///   call.
///
/// As an example, the following function might perform the position update for
/// a time-stepping code:
///
/// void update_position(T *data, const size_t count, const size_t offset,
///                     const E *env) {
///   for (size_t i = 0; i < count; ++i) {
///     data[i].position += data[i].velocity * env->delta_t;
///   }
/// }
///
/// For each action to be applied to a given sort of array, the user will need
/// to create another instance of an ArrayMapAction. These objects manage
/// registration of the action with the HPX-5 runtime, which means these objects
/// must be created before the call to dashmm::init().
template <typename T, typename E>
class ArrayMapAction {
 public:
  /// The function type for functions mapped onto Array elements.
  typedef void (*map_function_t)(T *, const size_t, const size_t, const E *);

  /// Construct the ArrayMapAction
  ///
  /// The only method for the object is the constructor, which takes a pointer
  /// to the function that will serve as the action represented by this object.
  /// This will ultimately register the needed information with the HPX-5
  /// runtime.
  ArrayMapAction(map_function_t f) {
    // If this is the first ArrayMapAction of this type, we need to register
    // the shared actions.
    if (!registered_) {
      HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_DEFAULT, root_, root_handler,
                          HPX_ACTION_T, HPX_POINTER, HPX_ADDR);
      HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_DEFAULT, spawn_, spawn_handler,
                          HPX_SIZE_T, HPX_SIZE_T, HPX_SIZE_T, HPX_SIZE_T,
                          HPX_ADDR, HPX_ACTION_T, HPX_POINTER, HPX_ADDR);
      registered_ = 1;
    }

    // Now register the specific function
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, leaf_, f,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER);
  }

 private:
  friend class Array<T>;

  static int root_handler(hpx_action_t leaf, const E *env,
                          hpx_addr_t meta_data) {
    ArrayMetaData *meta{nullptr};
    assert(hpx_gas_try_pin(meta_data, (void **)&meta));
    hpx_addr_t data = meta->data;
    size_t total_count = meta->count;
    hpx_gas_unpin(meta_data);

    // NOTE: Initially, for maximal "easiness" the library should pick the
    // decomposition factor. In the future, this should be extended and refined
    // in some way, where the user might give hints, or explicit instructions.
    size_t over_factor = 4;
    if (hpx_get_num_threads() == 1) {
      // Avoid pointless overhead for 1 thread
      over_factor = 1;
    }
    size_t n_per_chunk = total_count / (hpx_get_num_threads() * over_factor);
    size_t n_chunks{0};

    if (n_per_chunk == 0) {
      n_per_chunk = 1;
      n_chunks = total_count;
    } else {
      // Guarantees that all chunks, except possibly the last, have the same
      // size.
      n_chunks = total_count / n_per_chunk;
      if (total_count % n_per_chunk) {
        ++n_chunks;
      }
    }

    // Completion detection LCO
    hpx_addr_t alldone = hpx_lco_and_new(n_chunks);

    size_t offset = 0;
    hpx_call(HPX_HERE, spawn_, HPX_NULL, &total_count, &total_count,
             &offset, &n_per_chunk, &alldone, &leaf, &env, &data);

    hpx_lco_wait(alldone);
    hpx_lco_delete_sync(alldone);

    hpx_exit(HPX_SUCCESS);
  }

  static int spawn_handler(size_t count, size_t total_count, size_t offset,
                           size_t chunk_size, hpx_addr_t alldone,
                           hpx_action_t leaf, const E *env, hpx_addr_t data) {
    if (count <= chunk_size) {
      map_function_t lfunc = (map_function_t)hpx_action_get_handler(leaf);
      T *local{nullptr};
      assert(hpx_gas_try_pin(data, (void **)&local));
      lfunc(local, count, offset, env);
      hpx_gas_unpin(data);
      hpx_lco_set_lsync(alldone, 0, nullptr, HPX_NULL);
    } else {
      size_t num_chunks = count / chunk_size;
      if (count % chunk_size) {
        ++num_chunks;
      }

      size_t num_left = num_chunks / 2;
      size_t count_left = num_left * chunk_size;
      size_t count_right = count - count_left;
      size_t offset_right = offset + count_left;
      hpx_addr_t data_right = hpx_addr_add(data, sizeof(T) * count_left,
                                           sizeof(T) * total_count);

      hpx_call(HPX_HERE, spawn_, HPX_NULL, &count_left, &total_count, &offset,
               &chunk_size, &alldone, &leaf, &env, &data);
      hpx_call(HPX_HERE, spawn_, HPX_NULL, &count_right, &total_count,
               &offset_right, &chunk_size, &alldone, &leaf, &env, &data_right);
    }

    return HPX_SUCCESS;
  }

  static int registered_;
  static hpx_action_t root_;
  static hpx_action_t spawn_;
  hpx_action_t leaf_;
};

template <typename T, typename E>
int ArrayMapAction<T, E>::registered_ = 0;

template <typename T, typename E>
hpx_action_t ArrayMapAction<T, E>::root_ = HPX_ACTION_NULL;

template <typename T, typename E>
hpx_action_t ArrayMapAction<T, E>::spawn_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_ARRAY_MAP_ACTION_H__

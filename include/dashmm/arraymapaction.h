// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_ARRAY_MAP_ACTION_H__
#define __DASHMM_ARRAY_MAP_ACTION_H__


/// \file
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
/// a DASHMM array. It is a template requiring the type for the array
/// records, an environment that is passed to each invocation of the
/// associated action and an integer specifying how parallel to make the mapped
/// action.
///
/// To specify the action, the user will provide a function pointer with the
/// following signature: void func(T *, const size_t, const E *),
/// which has been aliased as map_function_t. When writing such a function,
/// it is important to know that the meaning of the provided arguments:
///
/// * The first is a pointer to an array of T objects. These are the data on
///   which the action will act.
/// * The second is the total count of records to be examined in the current
///   call of the function.
/// * The final is the environment passed into the Array<T>::map() function
///   call.
///
/// As an example, the following function might perform the position update for
/// a time-stepping code:
///
/// ~~~{.cc}
/// void update_position(T *data, const size_t count, const E *env) {
///   for (size_t i = 0; i < count; ++i) {
///     data[i].position += data[i].velocity * env->delta_t;
///   }
/// }
/// ~~~
///
/// For each action to be applied to a given sort of array, the user will need
/// to create another instance of an ArrayMapAction. These objects manage
/// registration of the action with the HPX-5 runtime, which means these objects
/// must be created before the call to dashmm::init().
///
/// The final template parameter gives the degree of parallelism to employ in
/// the resulting mapping of the action to the array. If the parameter is given
/// the value zero, a single thread per rank will process the entire array.
/// If the parameter is positive, then the array segment will be split into a
/// number of chunks equal to the number of HPX-5 scheduler threads times this
/// factor. Certain situations will decrease the number of chunks: if the
/// system is running with 1 thread, it will use only one chunk to avoid
/// pointless overhead; if there are too few records, a number of chunks will
/// be created equal to the number of records. This parameter has a default
/// value of 1; unless the amount of work per record is very non-uniform,
/// there should be little need to use another value.
template <typename T, typename E, int factor = 1>
class ArrayMapAction {
 public:
  static_assert(factor >= 0, "Decomposition factor must be non-negative");

  /// The function type for functions mapped onto Array elements.
  using map_function_t = void (*)(T *, const size_t, const E *);

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
                          HPX_SIZE_T, HPX_SIZE_T, HPX_SIZE_T,
                          HPX_ADDR, HPX_ACTION_T, HPX_POINTER, HPX_POINTER);
      registered_ = 1;
    }

    // Now register the specific function
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, leaf_, f,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER);
  }

 private:
  friend class Array<T>;

  /// Action that is the root of the parallel work spawn
  ///
  /// This action will start the parallel work. This takes care of a few
  /// tasks related to the management of the parallel spawn of the work.
  ///
  /// \param leaf - the action to take at the leaf computation
  /// \param env - the user supplied environment
  /// \param meta_data - global address of the Array's meta data
  ///
  /// \returns HPX_SUCCESS
  static int root_handler(hpx_action_t leaf, const E *env,
                          hpx_addr_t meta_data) {
    hpx_addr_t global = hpx_addr_add(meta_data,
                                     sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                     sizeof(ArrayMetaData));
    ArrayMetaData *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));
    char *data = local->data;
    size_t total_count = local->local_count;
    hpx_gas_unpin(global);

    if (total_count == 0) {
      hpx_exit(0, nullptr);
    }

    size_t over_factor = factor;
    if (factor == 0) {
      // Use only one chunk
      over_factor = 1;
    } else if (hpx_get_num_threads() == 1) {
      // Avoid pointless overhead for 1 thread
      over_factor = 1;
    } else {
      over_factor = hpx_get_num_threads() * factor;
    }
    size_t n_per_chunk = total_count / over_factor;
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

    hpx_call(HPX_HERE, spawn_, HPX_NULL, &total_count, &total_count,
             &n_per_chunk, &alldone, &leaf, &env, &data);

    hpx_lco_wait(alldone);
    hpx_lco_delete_sync(alldone);

    hpx_exit(0, nullptr);
  }

  /// Action serving the parallel spawn of the various chunks of work
  ///
  /// This will examine the inputs and spawn two more pieces of work if there
  /// are enough chunks assigned to this action.  Eventually, only one
  /// chunk will remain, and this will then initiate the leaf action on the
  /// chunk of records.
  ///
  /// \param count - the number of records for this action to map over
  /// \param total_count - the total number of records being mapped over
  /// \param chunk_size - the size in records of each chunk
  /// \param alldone - an LCO to manage completion detection
  /// \param leaf - the leaf action to take
  /// \param env - the environment to pass to the leaf action
  /// \param data - the global data to be handled by tis action
  ///
  /// \returns - HPX_SUCCESS
  static int spawn_handler(size_t count, size_t total_count,
                           size_t chunk_size, hpx_addr_t alldone,
                           hpx_action_t leaf, const E *env, char *data) {
    if (count <= chunk_size) {
      map_function_t lfunc = (map_function_t)hpx_action_get_handler(leaf);
      lfunc((T *)data, count, env);
      hpx_lco_set_lsync(alldone, 0, nullptr, HPX_NULL);
    } else {
      size_t num_chunks = count / chunk_size;
      if (count % chunk_size) {
        ++num_chunks;
      }

      size_t num_left = num_chunks / 2;
      size_t count_left = num_left * chunk_size;
      size_t count_right = count - count_left;
      char *data_right = data + sizeof(T) * count_left;

      hpx_call(HPX_HERE, spawn_, HPX_NULL, &count_left, &total_count,
               &chunk_size, &alldone, &leaf, &env, &data);
      hpx_call(HPX_HERE, spawn_, HPX_NULL, &count_right, &total_count,
               &chunk_size, &alldone, &leaf, &env, &data_right);
    }

    return HPX_SUCCESS;
  }

  static int registered_;
  static hpx_action_t root_;
  static hpx_action_t spawn_;
  hpx_action_t leaf_;
};

template <typename T, typename E, int factor>
int ArrayMapAction<T, E, factor>::registered_ = 0;

template <typename T, typename E, int factor>
hpx_action_t ArrayMapAction<T, E, factor>::root_ = HPX_ACTION_NULL;

template <typename T, typename E, int factor>
hpx_action_t ArrayMapAction<T, E, factor>::spawn_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_ARRAY_MAP_ACTION_H__

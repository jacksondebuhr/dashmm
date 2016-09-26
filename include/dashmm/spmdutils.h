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


#ifndef __DASHMM_BROADCAST_H__
#define __DASHMM_BROADCAST_H__


/// \file
/// \brief Ease of use class to provide value broadcast

#include "hpx/hpx.h"


namespace dashmm {


/// The broadcast value action identifier
extern hpx_action_t broadcast_value_action;


/// Broadcasts a given value from rank 0 to the other ranks
///
/// The provided value will be broadcast to each other rank. The provided
/// address will both supply the argument and receive the argument for rank
/// 0, and will receive the value for all other ranks.
///
/// The only requirement on the template parameter is that T needs to be
/// trivially copyable.
///
/// \param value - for rank 0, this is the input value; for all ranks, this
///                address will store the resulting value.
template <typename T>
void broadcast(T *value) {
  if (hpx_get_num_ranks() > 1) {
    hpx_run(&broadcast_value_action, value, value, sizeof(T));
  }
}


/// Get the rank of the calling locality
inline int get_my_rank() {
  return hpx_get_my_rank();
}


/// Get the number of ranks
inline int get_num_ranks() {
  return hpx_get_num_ranks();
}


} // dashmm


#endif // __DASHMM_BROADCAST_H__

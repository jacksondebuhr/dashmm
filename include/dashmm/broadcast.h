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


#include "hpx/hpx.h"


namespace dashmm {


extern hpx_action_t broadcast_value_action;


/// This broadcasts the given value from rank 0 to the other ranks
template <typename T>
void broadcast(T *value) {
  hpx_run(&broadcast_value_action, value, value, sizeof(T));
}


} // dashmm


#endif // __DASHMM_BROADCAST_H__
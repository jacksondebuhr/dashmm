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


#ifndef __DASHMM_REDUCTION_OPS_H__
#define __DASHMM_REDUCTION_OPS_H__


/// \file
/// \brief Action identifiers for common reduction operations


#include <hpx/hpx.h>


namespace dashmm {


/// Identity operation for integer summation
extern hpx_action_t int_sum_ident_op;

/// Operation for integer summation
extern hpx_action_t int_sum_op;

/// Identity operation for size_t summation
extern hpx_action_t size_sum_ident;

/// Operation for size_t summation
extern hpx_action_t size_sum_op;


} // namespace dashmm


#endif // __DASHMM_REDUCTION_OPS_H__

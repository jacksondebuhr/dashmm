// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
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

/// Identity operation for integer maximum
extern hpx_action_t int_max_ident_op;

/// Operation for integer maximum
extern hpx_action_t int_max_op;


} // namespace dashmm


#endif // __DASHMM_REDUCTION_OPS_H__

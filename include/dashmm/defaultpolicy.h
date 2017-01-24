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


#ifndef __DASHMM_DEFAULT_POLICY_H__
#define __DASHMM_DEFAULT_POLICY_H__


/// \file
/// \brief Specify default policies for use in the Library


#include "builtins/bhdistro.h"


namespace dashmm {


/// The default distribution policy to be used by Methods. This represents the
/// best all around policy that is available in DASHMM.
using DefaultDistributionPolicy = BHDistro;


} // namespace dashmm


#endif // __DASHMM_DEFAULT_POLICY_H__

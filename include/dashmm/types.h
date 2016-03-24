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


#ifndef __DASHMM_TYPES_H__
#define __DASHMM_TYPES_H__


/// \file include/dashmm/types.h
/// \brief Basic type definitions used throughout DASHMM


#include <complex>

#include "dashmm/domaingeometry.h"
#include "dashmm/index.h"
#include "dashmm/point.h"


namespace dashmm {


/// Return codes from DASHMM library calls
///
/// The possible return codes from DASHMM calls. For details about the
/// meanings of these, please see the individual library calls that generate
/// them.
enum ReturnCode {
  kSuccess = 0,
  kRuntimeError = 1,
  kIncompatible = 2,
  kAllocationError = 3,
  kInitError = 4,
  kFiniError = 5,
  kDomainError = 6
};


using dcomplex_t = std::complex<double>;


} // namespace dashmm


#endif // __DASHMM_TYPES_H__

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


#ifndef __DASHMM_TYPES_H__
#define __DASHMM_TYPES_H__


/// \file
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


/// Alias to standard library complex number type
using dcomplex_t = std::complex<double>;


/// Role codes for Expansion objects
///
/// Expansions are objects that can take up to four forms based on where they
/// are in the DAG. DAG nodes are associated with the nodes of either the source
/// or target tree. There are two roles that an Expansion can take for each tree
/// that the DAG node might be associated with: a primary expansion and an
/// intermediate expansion. For example, classic FMM uses only primary
/// expansion, which are called multipole expansions on the source tree, and
/// local expansions on the target tree. More advanced FMM implementations use
/// intermediate expansions on both trees. A fifth role is added to account for
/// cases where access is needed to the S->T operation, but no other expansion
/// data is required.
enum ExpansionRole {
  kSourcePrimary = 0,
  kSourceIntermediate = 1,
  kTargetPrimary = 2,
  kTargetIntermediate = 3,
  kNoRoleNeeded = 4
};

/// Operation codes to indicate the type of edge
enum class Operation {
  Nop,
  StoM,
  StoL,
  MtoM,
  MtoL,
  LtoL,
  MtoT,
  LtoT,
  StoT,
  MtoI,
  ItoI,
  ItoL
};


} // namespace dashmm


#endif // __DASHMM_TYPES_H__

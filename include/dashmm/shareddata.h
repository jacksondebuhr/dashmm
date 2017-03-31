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


#ifndef __DASHMM_SHAREDDATA_H__
#define __DASHMM_SHAREDDATA_H__


/// \file
/// \brief Interface to data shared across localities


#include "dashmm/domaingeometry.h"


namespace dashmm {

namespace shared {


/// Get the shared geometry data
///
/// \returns - Address of shared data
const DomainGeometry &geo();

/// Set the shared geometry data
///
/// \param g - the value to which to set the data
void set_geo(const DomainGeometry &g);


} // dashmm::shared

} // dashmm


#endif // __DASHMM_SHAREDDATA_H__
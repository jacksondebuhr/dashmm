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


#include "dashmm/shareddata.h"


namespace dashmm {

namespace shared {


/// The shared geometry data; there will be one per locality
DomainGeometry _shared_geo;


const DomainGeometry &geo() {
  return _shared_geo;
}


void set_geo(const DomainGeometry &g) {
  _shared_geo = g;
}


} // dashmm::shared


} // dashmm
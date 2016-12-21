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


#include "dashmm/shareddata.h"


namespace dashmm {

namespace shared {


/// The shared geometry data; there will be one per locality
DomainGeometry _shared_geo;


/*
NOTE: At the moment, we rely on having written all the use cases to be
sure that this is a safe thing to use. In the future, we likely want to make
this more resilient to multiple threads, and timing issues. This would
potentially mean an LCO to protect it, and a flag to indicate the status of
this data.

TODO: Consider all of that.
*/


const DomainGeometry &geo() {
  return _shared_geo;
}


void set_geo(const DomainGeometry &g) {
  _shared_geo = g;
}


} // dashmm::shared


} // dashmm
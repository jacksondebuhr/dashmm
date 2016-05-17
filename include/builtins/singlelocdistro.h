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


#ifndef __DASHMM_SINGLE_LOC_DISTRO_H__
#define __DASHMM_SINGLE_LOC_DISTRO_H__


namespace dashmm {


// A distribution policy that will place everything on one locality
// only use if it definitely fits in memory
class SingleLocality {
public:
  void compute_distribution();
}


} // dashmm


#endif // __DASHMM_SINGLE_LOC_DISTRO_H__

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


#include "dashmm/dag.h"


namespace dashmm {


/// This distribution policy places all content on the root locality
///
/// This is intended to be the simplest possibly distribution policy. It is
/// unlikely to be a good choice, unless one is using only one locality to
/// begin with.
class SingleLocality {
 public:
  SingleLocality(int loc = 0) : locality_{loc} { }

  void compute_distribution(DAG &dag);

 private:
  int locality_;
};


} // dashmm


#endif // __DASHMM_SINGLE_LOC_DISTRO_H__

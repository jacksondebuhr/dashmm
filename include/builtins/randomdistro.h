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


#ifndef __DASHMM_RANDOM_DISTRO_H__
#define __DASHMM_RANDOM_DISTRO_H__


#include "dashmm/dag.h"


namespace dashmm {


class RandomDistro {
 public:
  RandomDistro(int seed = 137) : seed_{seed} { }
  void compute_distribution(DAG &dag);

 private:
  int seed_;
};


} // dashmm


#endif // __DASHMM_RANDOM_DISTRO_H__

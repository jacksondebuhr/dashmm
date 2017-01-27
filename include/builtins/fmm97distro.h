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


#ifndef __DASHMM_FMM97_DISTRO_H__
#define __DASHMM_FMM97_DISTRO_H__


/// \file
/// \brief definition of FMM97 Distribution Policy


#include "dashmm/dag.h"


namespace dashmm {


class FMM97Distro {
 public:
  FMM97Distro() { }
  void compute_distribution(DAG &dag);
  static void assign_for_source(DAGInfo &dag, int locality);
  static void assign_for_target(DAGInfo &dag, int locality);

 private:
  static bool color_comparison(const DAGNode *a, const DAGNode *b);
  static void compute_locality(DAGNode *node);
};


} // dashmm


#endif // __DASHMM_FMM97_DISTRO_H__

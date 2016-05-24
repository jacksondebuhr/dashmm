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


#include "builtins/singlelocdistro.h"


namespace dashmm {


void SingleLocality::compute_distribution(
    const SharedData<DomainGeometry> &domain,
    const std::vector<DAGNode *> &sources,
    const std::vector<DAGNode *> &targets,
    const std::vector<DAGNode *> &internal) {
  auto b = internal.begin();
  auto e = internal.end();
  for (auto i = b; i != e; ++i) {
    (*i)->locality = 0;
  }
}


}

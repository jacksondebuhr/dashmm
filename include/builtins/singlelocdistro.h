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


#include <vector>

#include "dashmm/daginfo.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/shareddata.h"


namespace dashmm {


/// This distribution policy places all content on the root locality
///
/// This is intended to be the simplest possibly distribution policy. It is
/// unlikely to be a good choice, unless one is using only one locality to
/// begin with.
class SingleLocality {
public:
  void compute_distribution(const SharedData<DomainGeometry> &domain,
                            const std::vector<DAGNode *> &sources,
                            const std::vector<DAGNode *> &targets,
                            const std::vector<DAGNode *> &internal);
};


} // dashmm


#endif // __DASHMM_SINGLE_LOC_DISTRO_H__

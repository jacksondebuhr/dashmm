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


/// \file
/// \brief Implementation of JSON format DAG output
///
/// The intent of this file it to make an easily digestible form of the
/// DAG information for use in visualization tools. This implements JSON
/// formatted data, as it is human readable (if boring) and because there
/// are quality JSON readers in a wide array of languages and frameworks.
///
/// NOTE: This is not as robust as other portions of DASHMM, and should be
/// considered to be experimental. The default mode of operation of DASHMM
/// will not even call into these routines, and must be manually enabled
/// by modifying the source code of the library. Please note the intentional
/// vagueness.


#include "dashmm/dag.h"


namespace dashmm {


Index DAGNode::index() const {
  return parent_->index();
}

bool DAGNode::is_parts() const {
  return parent_->parts() == this;
}

bool DAGNode::is_normal() const {
  return parent_->normal() == this;
}

bool DAGNode::is_interm() const {
  return parent_->interm() == this;
}


} // dashmm

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


#ifndef __DASHMM_EVALUATOR_H__
#define __DASHMM_EVALUATOR_H__


template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename, typename> class Method>
class Evaluator {
 public:
  Evaluator() {
    // TODO: register all actions in the system
  }
 private:
};


#endif // __DASHMM_EVALUATOR_H__

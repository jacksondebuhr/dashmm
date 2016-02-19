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


#ifndef __DASHMM_DIRECT_METHOD_H__
#define __DASHMM_DIRECT_METHOD_H__


/// \file include/direct_method.h
/// \brief Direction summation method


#include "include/expansionlco.h"
#include "include/method.h"
#include "include/node.h"


namespace dashmm {


template <typename Source, typename Target,
          template <typename, typename> class Expansion>
class Direct {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Direct<Source, Target, Expansion>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Direct>;
  using sourcenode_t = SourceNode<Source, Target, Expansion, Direct>;
  using targetnode_t = TargetNode<Source, Target, Expansion, Direct>;

  void generate(sourcenode_t &curr, int n_digits) const {
    curr.set_expansion(std::unique_ptr<expansion_t>{
        new expansion_t{Point{0.0, 0.0, 0.0}, n_digits}
      });
  }

  void aggregate(sourcenode_t &curr, int n_digits) const { }

  void inherit(targetnode_t &curr, int n_digits, size_t which_child) const {
    curr.set_expansion(std::unique_ptr<expansion_t>{
        new expansion_t{Point{0.0, 0.0, 0.0}, n_digits}
      });
  }

  void process(targetnode_t &curr, std::vector<sourcenode_t> &consider,
               bool curr_is_leaf) const {
    std::vector<sourcenode_t> newcons{ };
    do {
      for (auto i = consider.begin(); i != consider.end(); ++i) {
        if (i->is_leaf()) {
          if (curr_is_leaf) {
            expansionlco_t expand = i->expansion();
            targetlco_t targets = curr.parts();
            sourceref_t sources = i->parts();
            expand.S_to_T(sources, targets);
          } else {
            newcons.push_back(*i);
          }
        } else {
          for (size_t j = 0; j < 8; ++j) {
            sourcenode_t kid = i->child(j);
            if (kid.is_valid()) {
              newcons.push_back(kid);
            }
          }
        }
      }

      consider = std::move(newcons);
      newcons.clear();
    } while (curr_is_leaf && !consider.empty());
  }

  // One option would be to cause this test to always return false, meaning the
  // trees would not refine at all, and thus create two huge nodes, and then
  // run S->T for the given expansion. However, we can get some parallelism if
  // we go ahead and refine where we can. Just because we are using a dumb
  // method, does not mean we have to do that dumb thing in serial.
  bool refine_test(bool same_sources_and_targets, const targetnode_t &curr,
                   const std::vector<sourcenode_t> &consider) const {
    return true;
  }
};


} // namespace dashmm


#endif // __DASHMM_DIRECT_METHOD_H__

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


#ifndef __DASHMM_DIRECT_METHOD_H__
#define __DASHMM_DIRECT_METHOD_H__


/// \file
/// \brief Direction summation method


#include <cstdlib>

#include <vector>

#include "dashmm/arrayref.h"
#include "dashmm/defaultpolicy.h"
#include "dashmm/expansionlco.h"
#include "dashmm/targetlco.h"
#include "dashmm/tree.h"

#include "builtins/bhdistro.h"


namespace dashmm {


/// A Method to implement the direct summation method
///
/// This is an O(N^2) method, so it should be used for only small numbers
/// of targets.
template <typename Source, typename Target,
          template <typename, typename> class Expansion>
class Direct {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Direct<Source, Target, Expansion>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;

  using distropolicy_t = BHDistro;

  Direct() { }

  /// In generate, Direct does nothing.
  void generate(sourcenode_t *curr, DomainGeometry *domain) const {
    curr->dag.add_parts();
  }

  /// In aggregate, Direct does nothing.
  void aggregate(sourcenode_t *curr, DomainGeometry *domain) const { }

  /// In inherit, Direct does nothing.
  void inherit(targetnode_t *curr, DomainGeometry *domain,
               bool curr_is_leaf) const {
    if (curr_is_leaf) {
      curr->dag.add_parts();
    }
  }

  /// In process, Direct will collect all leaf nodes and apply S->T for each
  void process(targetnode_t *curr, std::vector<sourcenode_t *> &consider,
               bool curr_is_leaf, DomainGeometry *domain) const {
    std::vector<sourcenode_t *> newcons{ };
    do {
      for (auto i = consider.begin(); i != consider.end(); ++i) {
        if ((*i)->is_leaf()) {
          if (curr_is_leaf) {
            curr->dag.StoT(&(*i)->dag,
                           expansion_t::weight_estimate(Operation::StoT));
          } else {
            newcons.push_back(*i);
          }
        } else {
          for (size_t j = 0; j < 8; ++j) {
            sourcenode_t *kid = (*i)->child[j];
            if (kid != nullptr) {
              newcons.push_back(kid);
            }
          }
        }
      }

      consider = std::move(newcons);
      newcons.clear();
    } while (curr_is_leaf && !consider.empty());
  }

  /// Refinement test for Direct method.
  ///
  /// One option would be to cause this test to always return false, meaning the
  /// trees would not refine at all, and thus create two huge nodes, and then
  /// run S->T for the given expansion. However, we can get some parallelism if
  /// we go ahead and refine where we can. Just because we are using a dumb
  /// method, does not mean we have to do that dumb thing in serial.
  bool refine_test(bool same_sources_and_targets, const targetnode_t *curr,
                   const std::vector<sourcenode_t *> &consider) const {
    return true;
  }
};


} // namespace dashmm


#endif // __DASHMM_DIRECT_METHOD_H__

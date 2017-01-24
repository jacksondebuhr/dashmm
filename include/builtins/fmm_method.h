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


#ifndef __DASHMM_FMM_METHOD_H__
#define __DASHMM_FMM_METHOD_H__


/// \file
/// \brief Declaration of FMM Method


#include "dashmm/arrayref.h"
#include "dashmm/defaultpolicy.h"
#include "dashmm/expansionlco.h"
#include "dashmm/index.h"
#include "dashmm/targetlco.h"
#include "dashmm/tree.h"

#include "builtins/bhdistro.h"


namespace dashmm {


/// A Method to implement classic FMM
///
template <typename Source, typename Target,
          template <typename, typename> class Expansion>
class FMM {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = FMM<Source, Target, Expansion>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;

  using distropolicy_t = BHDistro;

  FMM() { }

  void generate(sourcenode_t *curr, DomainGeometry *domain) const {
    curr->dag.add_parts();
    if (curr->idx.level() >= 2) {
      // If \p curr is of level 0 or 1, \p curr is not well separated
      // from any node. As a result, there is no need to create the
      // normal expansion.
      assert(curr->dag.add_normal() == true);
      curr->dag.StoM(&curr->dag,
                     expansion_t::weight_estimate(Operation::StoM));
    }
  }

  void aggregate(sourcenode_t *curr, DomainGeometry *domain) const {
    if (curr->idx.level() >= 2) {
      assert(curr->dag.add_normal() == true);
      for (size_t i = 0; i < 8; ++i) {
        sourcenode_t *kid = curr->child[i];
        if (kid != nullptr) {
          curr->dag.MtoM(&kid->dag,
                         expansion_t::weight_estimate(Operation::MtoM));
        }
      }
    }
  }

  void inherit(targetnode_t *curr, DomainGeometry *domain,
               bool curr_is_leaf) const {
    if (curr_is_leaf) {
      curr->dag.add_parts();
    }

    if (curr->idx.level() >= 2) {
      assert(curr->dag.add_normal() == true);

      if (curr->idx.level() >= 3) {
        curr->dag.LtoL(&curr->parent->dag,
                       expansion_t::weight_estimate(Operation::LtoL));
      }
    }
  }

  void process(targetnode_t *curr, std::vector<sourcenode_t *> &consider,
               bool curr_is_leaf, DomainGeometry *domain) const {
    Index t_index = curr->idx;
    int StoL = expansion_t::weight_estimate(Operation::StoL);
    int StoT = expansion_t::weight_estimate(Operation::StoT);
    int LtoT = expansion_t::weight_estimate(Operation::LtoT);
    int MtoL = expansion_t::weight_estimate(Operation::MtoL);

    if (curr_is_leaf) {
      for (auto S = consider.begin(); S != consider.end(); ++S) {
        if ((*S)->idx.level() < t_index.level()) {
          if (well_sep_test_asymmetric(t_index, (*S)->idx)) {
            curr->dag.StoL(&(*S)->dag, StoL);
          } else {
            curr->dag.StoT(&(*S)->dag, StoT);
          }
        } else {
          if (well_sep_test((*S)->idx, t_index)) {
            curr->dag.MtoL(&(*S)->dag, MtoL);
          } else {
            proc_coll_recur(curr, *S);
          }
        }
      }

      if (curr->idx.level() >= 2)
        curr->dag.LtoT(&curr->dag, LtoT);
    } else {
      std::vector<sourcenode_t *> newcons{ };

      for (auto S = consider.begin(); S != consider.end(); ++S) {
        if ((*S)->idx.level() < t_index.level()) {
          if (well_sep_test_asymmetric(t_index, (*S)->idx)) {
            curr->dag.StoL(&(*S)->dag, StoL);
          } else {
            newcons.push_back(*S);
          }
        } else {
          if (well_sep_test((*S)->idx, t_index)) {
            curr->dag.MtoL(&(*S)->dag, MtoL);
          } else {
            bool S_is_leaf = true;
            for (size_t i = 0; i < 8; ++i) {
              sourcenode_t *child = (*S)->child[i];
              if (child != nullptr) {
                newcons.push_back(child);
                S_is_leaf = false;
              }
            }

            if (S_is_leaf) {
              newcons.push_back(*S);
            }
          }
        }
      }

      consider = std::move(newcons);
    }
  }

  bool refine_test(bool same_sources_and_targets, const targetnode_t *curr,
                   const std::vector<sourcenode_t *> &consider) const {
    if (same_sources_and_targets) {
      return true;
    }

    for (auto i = consider.begin(); i != consider.end(); ++i) {
      if ((*i)->idx.level() == curr->idx.level()) {
        if (!well_sep_test((*i)->idx, curr->idx) && !(*i)->is_leaf()) {
          return true;
        }
      }
    }

    return false;
  }

  bool well_sep_test_asymmetric(Index smaller, Index larger) const {
    int delta = smaller.level() - larger.level();
    int shift = (1 << delta) - 1;
    Index l_mod{larger.x() << delta, larger.y() << delta,
          larger.z() << delta, smaller.level()};
    Index l_mod_top{l_mod.x() + shift, l_mod.y() + shift, l_mod.z() + shift,
          l_mod.level()};

    // NOTE: shortcut evaluation
    return (l_mod.x() - smaller.x() > 1) || (smaller.x() - l_mod_top.x() > 1)
        || (l_mod.y() - smaller.y() > 1) || (smaller.y() - l_mod_top.y() > 1)
        || (l_mod.z() - smaller.z() > 1) || (smaller.z() - l_mod_top.z() > 1);
  }

  bool well_sep_test(Index source, Index target) const {
    // When the nodes are the same level, we just need to have at least one
    // index that is different by +/- 2 or more.
    if (abs(source.x() - target.x()) > 1) return true;
    if (abs(source.y() - target.y()) > 1) return true;
    if (abs(source.z() - target.z()) > 1) return true;
    return false;
  }

  void proc_coll_recur(targetnode_t *T, sourcenode_t *S) const {
    int MtoT = expansion_t::weight_estimate(Operation::MtoT);
    int StoT = expansion_t::weight_estimate(Operation::StoT);

    if (well_sep_test_asymmetric(S->idx, T->idx)) {
      S->dag.MtoT(&T->dag, MtoT);
    } else {
      if (S->is_leaf()) {
        T->dag.StoT(&S->dag, StoT);
      } else {
        for (size_t i = 0; i < 8; ++i) {
          sourcenode_t *child = S->child[i];
          if (child != nullptr) {
            proc_coll_recur(T, child);
          }
        }
      }
    }
  }
};


} // namespace dashmm


#endif // __DASHMM_FMM_METHOD_H__

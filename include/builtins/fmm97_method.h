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

#ifndef __DASHMM_FMM97_METHOD_H__
#define __DASHMM_FMM97_METHOD_H__

/// \file include/builtins/fmm97_method.h
/// \brief Declaration of FMM method with merge and shift 

#include "dashmm/arrayref.h"
#include "dashmm/defaultpolicy.h"
#include "dashmm/expansionlco.h"
#include "dashmm/index.h"
#include "dashmm/targetlco.h"
#include "dashmm/tree.h"

namespace dashmm {

/// A Method to implement FMM with merge and shift
///
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          typename DistroPolicy = DefaultDistributionPolicy>
class FMM97 {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = FMM97<Source, Target, Expansion, DistroPolicy>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, FMM97,
                                      DistroPolicy>;
  using sourceref_t = ArrayRef<Source>;
  using sourcenode_t = TreeNode<Source, Target, Source, Expansion, FMM97,
                                DistroPolicy>;
  using targetnode_t = TreeNode<Source, Target, Target, Expansion, FMM97,
                                DistroPolicy>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, FMM97, DistroPolicy>;

  FMM97() {} 

  void generate(sourcenode_t *curr, DomainGeometry *domain) const {
    curr->dag.add_parts();
    if (curr->idx.level() >= 2) {
      // If \p curr is of level 0 or 1, \p curr is not well separated
      // from any target node. As a result, there is no need to create
      // the normal or intermediate expansion 
      curr->dag.add_normal();
      curr->dag.add_interm(); 
      curr->dag.StoM(&curr->dag);
      curr->dag.MtoI(&curr->dag); 
    }
  }

  void aggregate(sourcenode_t *curr, DomainGeometry *domain) const {
    if (curr->idx.level() >= 2) {
      curr->dag.add_normal();
      curr->dag.add_interm(); 
      curr->dag.MtoI(&curr->dag); 
      for (size_t i = 0; i < 8; ++i) {
        sourcenode_t *kid = curr->child[i];
        if (kid != nullptr) {
          curr->dag.MtoM(&kid->dag);
        }
      }
    }
  }

  void inherit(targetnode_t *curr, DomainGeometry *domain,
               bool curr_is_leaf) const {
    if (curr_is_leaf) 
      curr->dag.add_parts(); 

    if (curr->idx.level() >= 2) {
      curr->dag.add_normal(); 
      curr->dag.ItoL(&curr->parent->dag); 
      
      if (curr->idx.level() >= 3) 
        curr->dag.LtoL(&curr->parent->dag); 
    }
  }

  void process(targetnode_t *curr, std::vector<sourcenode_t *> &consider,
               bool curr_is_leaf, DomainGeometry *domain) const {
    Index t_index = curr->idx;

    // If a source node S is at the same level as \p curr and is well separated
    // from \p curr, skip S as its contribution to \p curr has been processed by
    // the parent of \p curr. 
    if (curr_is_leaf) {
      for (auto S = consider.begin(); S != consider.end(); ++S) {
        if ((*S)->idx.level() < t_index.level()) {
          if (well_sep_test_asymmetric(t_index, (*S)->idx)) {
            curr->dag.StoL(&(*S)->dag);
          } else {
            curr->dag.StoT(&(*S)->dag);
          }
        } else {
          if (!well_sep_test((*S)->idx, t_index)) {
            proc_coll_recur(curr, *S);
          }
        }
      }

      if (curr->idx.level() >= 2) 
        curr->dag.LtoT(&curr->dag);
    } else {
      if (curr->idx.level() >= 1) 
        curr->dag.add_interm(); 

      std::vector<sourcenode_t *> newcons{};

      for (auto S = consider.begin(); S != consider.end(); ++S) {
        if ((*S)->idx.level() < t_index.level()) {
          if (well_sep_test_asymmetric(t_index, (*S)->idx)) {
            curr->dag.StoL(&(*S)->dag);
          } else {
            newcons.push_back(*S);
          }
        } else {
          if (!well_sep_test((*S)->idx, t_index)) {
            bool do_I2I = (curr->idx.level() >= 1 && 
                           ((*S)->idx.x() != t_index.x() || 
                            (*S)->idx.y() != t_index.y() ||
                            (*S)->idx.z() != t_index.z())); 

            bool S_is_leaf = true; 

            for (size_t i = 0; i < 8; ++i) {
              sourcenode_t *child = (*S)->child[i]; 
              if (child != nullptr) {
                newcons.push_back(child); 
                S_is_leaf = false; 
                if (do_I2I) 
                  curr->dag.ItoI(&child->dag); 
              }
            }

            if (S_is_leaf) 
              newcons.push_back(*S); 
          }
        }
      } // endfor

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
    if (well_sep_test_asymmetric(S->idx, T->idx)) {
      S->dag.MtoT(&T->dag);
    } else {
      if (S->is_leaf()) {
        T->dag.StoT(&S->dag);
      } else {
        for (size_t i = 0; i < 8; ++i) {
          sourcenode_t *child = S->child[i];
          if (child != nullptr)
            proc_coll_recur(T, child);
        }
      }
    }
  }
};


} // namespace dashmm


#endif // __DASHMM_FMM_METHOD_H__

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


#ifndef __DASHMM_FMM_METHOD_H__
#define __DASHMM_FMM_METHOD_H__


/// \file include/builtins/fmm_method.h
/// \brief Declaration of FMM Method


#include "dashmm/expansionlco.h"
#include "dashmm/index.h"
#include "dashmm/sourcenode.h"
#include "dashmm/sourceref.h"
#include "dashmm/targetlco.h"
#include "dashmm/targetnode.h"


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
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, FMM>;
  using sourceref_t = SourceRef<Source>;
  using sourcenode_t = SourceNode<Source, Target, Expansion, FMM>;
  using targetnode_t = TargetNode<Source, Target, Expansion, FMM>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, FMM>;

  void generate(sourcenode_t &curr, int n_digits) const {
    curr.set_expansion(std::unique_ptr<expansion_t>{
        new expansion_t{curr.center(), n_digits}
      });
    expansionlco_t currexp = curr.expansion();
    sourceref_t sources = curr.parts();
    double scale = 1.0 / curr.size();
    currexp.S_to_M(curr.center(), sources, scale);
  }

  void aggregate(sourcenode_t &curr, int n_digits) const {
    curr.set_expansion(std::unique_ptr<expansion_t>{
        new expansion_t{curr.center(), n_digits}
      });
    expansionlco_t currexp = curr.expansion();
    for (size_t i = 0; i < 8; ++i) {
      sourcenode_t kid = curr.child(i);
      if (kid.is_valid()) {
        expansionlco_t kexp = kid.expansion();
        currexp.M_to_M(kexp, i, kid.size());
      }
    }
  }

  void inherit(targetnode_t &curr, int n_digits, size_t which_child) const {
    curr.set_expansion(std::unique_ptr<expansion_t>{
        new expansion_t{curr.center(), n_digits}
      });
    expansionlco_t currexp = curr.expansion();

    if (curr.parent().is_valid()) {
      expansionlco_t pexp = curr.parent().expansion();
      currexp.L_to_L(pexp, which_child, curr.size());
    }
  }

  void process(targetnode_t &curr, std::vector<sourcenode_t> &consider,
               bool curr_is_leaf) const {
    expansionlco_t currexp = curr.expansion();
    double scale = 1.0 / curr.size();
    targetlco_t targets = curr.parts();
    Index t_index = curr.index();

    if (curr_is_leaf) {
      for (auto S = consider.begin(); S != consider.end(); ++S) {
        if (S->level() < curr.level()) {
          if (well_sep_test_asymmetric(t_index, S->index())) {
            sourceref_t sources = S->parts();
            currexp.S_to_L(curr.center(), sources, scale);
          } else {
            sourceref_t sources = S->parts();
            expansionlco_t expand = S->expansion();
            expand.S_to_T(sources, targets);
          }
        } else {
          if (well_sep_test(S->index(), curr.index())) {
            expansionlco_t expand = S->expansion();
            double s_size = S->size();
            Index s_index = S->index();
            currexp.M_to_L(expand, s_index, s_size, t_index);
          } else {
            proc_coll_recur(curr, *S);
          }
        }
      }

      currexp.L_to_T(targets, scale);
    } else {
      std::vector<sourcenode_t> newcons{ };

      for (auto S = consider.begin(); S != consider.end(); ++S) {
        if (S->level() < curr.level()) {
          if (well_sep_test_asymmetric(t_index, S->index())) {
            sourceref_t sources = S->parts();
            currexp.S_to_L(curr.center(), sources, scale);
          } else {
            // a 1; add it to newcons
            newcons.push_back(*S);
          }
        } else {
          if (well_sep_test(S->index(), curr.index())) {
            expansionlco_t expand = S->expansion();
            double s_size = S->size();
            Index s_index = S->index();
            currexp.M_to_L(expand, s_index, s_size, t_index);
          } else {
            bool S_is_leaf = true;
            for (size_t i = 0; i < 8; ++i) {
              sourcenode_t child = S->child(i);
              if (child.is_valid()) {
                newcons.push_back(child);
                S_is_leaf = false;
              }
            }

            if (S_is_leaf)
              newcons.push_back(*S);
          }
        }
      }

      consider = std::move(newcons);
    }
  }

  bool refine_test(bool same_sources_and_targets, const targetnode_t &curr,
                   const std::vector<sourcenode_t> &consider) const {
    if (same_sources_and_targets) {
      return true;
    }

    for (auto i = consider.begin(); i != consider.end(); ++i) {
      if (i->level() == curr.level()) {
        if (!well_sep_test(i->index(), curr.index()) && !i->is_leaf()) {
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

    // in the asymmetric case, we need to check both sides of the larger box
    // to make sure there is a gap of one of the smaller boxes
    if (abs(smaller.x() - l_mod.x()) > 1
        && abs(smaller.x() - l_mod_top.x()) > 1) return true;
    if (abs(smaller.y() - l_mod.y()) > 1
        && abs(smaller.y() - l_mod_top.y()) > 1) return true;
    if (abs(smaller.z() - l_mod.z()) > 1
        && abs(smaller.z() - l_mod_top.z()) > 1) return true;
    return false;
  }

  bool well_sep_test(Index source, Index target) const {
    // When the nodes are the same level, we just need to have at least one
    // index that is different by +/- 2 or more.
    if (abs(source.x() - target.x()) > 1) return true;
    if (abs(source.y() - target.y()) > 1) return true;
    if (abs(source.z() - target.z()) > 1) return true;
    return false;
  }

  void proc_coll_recur(targetnode_t &T, sourcenode_t &S) const {
    if (well_sep_test_asymmetric(S.index(), T.index())) {
      expansionlco_t expand = S.expansion();
      targetlco_t targets = T.parts();
      double scale = 1.0 / S.size();
      expand.M_to_T(targets, scale);
    } else {
      if (S.is_leaf()) {
        expansionlco_t expand = S.expansion();
        targetlco_t targets = T.parts();
        sourceref_t sources = S.parts();
        expand.S_to_T(sources, targets);
      } else {
        for (size_t i = 0; i < 8; ++i) {
          sourcenode_t child = S.child(i);
          if (child.is_valid())
            proc_coll_recur(T, child);
        }
      }
    }
  }
};


} // namespace dashmm


#endif // __DASHMM_FMM_METHOD_H__

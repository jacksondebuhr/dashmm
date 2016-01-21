/// \file src/fmm_method.cc
/// \brief Implementation of FMM Method

#include "include/fmm_method.h"

namespace dashmm {

void FMM::generate(SourceNode &curr, const ExpansionRef expand) const {
  int n_digits = expand.accuracy(); 
  curr.set_expansion(expand.get_new_expansion(curr.center(), n_digits));
  ExpansionRef currexp = curr.expansion(); 
  SourceRef sources = curr.parts(); 
  double scale = 1.0 / curr.size(); 
  currexp.S_to_M(curr.center(), sources, scale); 
}

void FMM::aggregate(SourceNode &curr, const ExpansionRef expand) const {
  int n_digits = expand.accuracy(); 
  curr.set_expansion(expand.get_new_expansion(curr.center(), n_digits)); 
  ExpansionRef currexp = curr.expansion(); 

  for (size_t i = 0; i < 8; ++i) {
    SourceNode kid = curr.child(i); 
    if (kid.is_valid()) {
      ExpansionRef kexp = kid.expansion(); 
      currexp.M_to_M(kexp, i, kid.size());
    }
  }
}

void FMM::inherit(TargetNode &curr, const ExpansionRef expand, 
                  size_t which_child) const {
  int n_digits = expand.accuracy(); 
  curr.set_expansion(expand.get_new_expansion(curr.center(), n_digits)); 
  ExpansionRef currexp = curr.expansion(); 

  if (curr.parent().is_valid()) {
    ExpansionRef pexp = curr.parent().expansion(); 
    currexp.L_to_L(pexp, which_child, curr.size()); 
  } else {
    curr.set_expansion(expand.get_new_expansion(curr.center(), n_digits));
  } 
}

void FMM::process(TargetNode &curr, std::vector<SourceNode> &consider, 
                  bool curr_is_leaf) const {
  ExpansionRef currexp = curr.expansion(); 
  double scale = 1.0 / curr.size(); 
  TargetRef targets = curr.parts(); 
  Index t_index = curr.index(); 

  if (curr_is_leaf) {
    for (auto S = consider.begin(); S != consider.end(); ++S) {
      if (S->level() < curr.level()) {
        if (well_sep_test_asymmetric(t_index, S->index())) {
          SourceRef sources = S->parts(); 
          currexp.S_to_L(curr.center(), sources, scale); 
        } else {
          SourceRef sources = S->parts(); 
          ExpansionRef expand = S->expansion(); 
          expand.S_to_T(sources, targets);
        }
      } else {
        if (well_sep_test(S->index(), curr.index())) {
          ExpansionRef expand = S->expansion();           
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
    std::vector<SourceNode> newcons{}; 
    
    for (auto S = consider.begin(); S != consider.end(); ++S) {
      if (S->level() < curr.level()) {
        if (well_sep_test_asymmetric(t_index, S->index())) {
          SourceRef sources = S->parts(); 
          currexp.S_to_L(curr.center(), sources, scale); 
        } else {
          // a 1; add it to newcons
          newcons.push_back(*S); 
        }
      } else {
        if (well_sep_test(S->index(), curr.index())) {
          ExpansionRef expand = S->expansion(); 
          double s_size = S->size(); 
          Index s_index = S->index(); 
          currexp.M_to_L(expand, s_index, s_size, t_index); 
        } else {
          bool S_is_leaf = true; 
          for (size_t i = 0; i < 8; ++i) {
            SourceNode child = S->child(i); 
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

bool FMM::refine_test(bool same_sources_and_targets, const TargetNode &curr,
                      const std::vector<SourceNode> &consider) const {
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

bool FMM::well_sep_test_asymmetric(Index smaller, Index larger) const {
  int delta = smaller.level() - larger.level();
  int shift = (1 << delta) - 1;
  Index l_mod{larger.x() << delta, larger.y() << delta,
              larger.z() << delta, smaller.level()};
  Index l_mod_top{l_mod.x() + shift, l_mod.y() + shift, l_mod.z() + shift,
                  l_mod.level()};

  //in the asymmetric case, we need to check both sides of the larger box
  // to make sure there is a gap of one of the smaller boxes
  if (abs(smaller.x() - l_mod.x()) > 1
       && abs(smaller.x() - l_mod_top.x()) > 1) return true;
  if (abs(smaller.y() - l_mod.y()) > 1
       && abs(smaller.y() - l_mod_top.y()) > 1) return true;
  if (abs(smaller.z() - l_mod.z()) > 1
       && abs(smaller.z() - l_mod_top.z()) > 1) return true;
  return false;
}

bool FMM::well_sep_test(Index source, Index target) const {
  //When the nodes are the same level, we just need to have at least one
  // index that is different by +/- 2 or more.
  if (abs(source.x() - target.x()) > 1) return true;
  if (abs(source.y() - target.y()) > 1) return true;
  if (abs(source.z() - target.z()) > 1) return true;
  return false;
}

void FMM::proc_coll_recur(TargetNode &T, SourceNode &S) const {
  if (well_sep_test_asymmetric(S.index(), T.index())) {
    ExpansionRef expand = S.expansion(); 
    TargetRef targets = T.parts(); 
    double scale = 1.0 / T.size(); 
    expand.M_to_T(targets, scale); 
  } else {
    if (S.is_leaf()) {
      ExpansionRef expand = S.expansion(); 
      TargetRef targets = T.parts(); 
      SourceRef sources = S.parts(); 
      expand.S_to_T(sources, targets); 
    } else {
      for (size_t i = 0; i < 8; ++i) {
        SourceNode child = S.child(i); 
        if (child.is_valid()) 
          proc_coll_recur(T, child);
      }
    }
  }
}


} // namespace dashmm

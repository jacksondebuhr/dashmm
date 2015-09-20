#include "FMM_method.h"

#include <cassert>

#include "index.h"


namespace dashmm {


void FMM_Method::generate(SourceNode *curr, const Expansion *expand) const {
  curr->set_expansion(expand->S_to_M(curr->center(), curr->first(),
                                     curr->last()));
}


void FMM_Method::aggregate(SourceNode *curr, const Expansion *expand) const {
  curr->set_expansion(expand->get_new_expansion(curr->center()));

  for (size_t i = 0; i < 8; ++i) {
    SourceNode *kid = curr->child(i);
    if (kid) {
      std::unique_ptr<Expansion> term{
                      kid->expansion()->M_to_M(i, kid->size())};
      curr->expansion()->add_expansion(term.get());
    }
  }
}


void FMM_Method::inherit(TargetNode *curr, const Expansion *proto,
                         size_t which_child) const {
  if (curr->parent() == nullptr) {
    curr->set_expansion(proto->get_new_expansion(curr->center()));
  } else {
    curr->set_expansion(
        curr->parent()->expansion()->L_to_L(which_child, curr->size()));
  }
}


void FMM_Method::process(TargetNode *curr, std::vector<SourceNode *> &consider,
                         bool curr_is_leaf) const {
  std::vector<SourceNode *> newcons{};

  for (auto source = consider.begin(); source != consider.end(); ++source) {
    //everything in consider should be curr's size or larger.
    assert((*source)->level() <= curr->level());

    if ((*source)->level() < curr->level()) {
      assert((*source)->is_leaf());
      //a 1 or a 4
      if (well_sep_test_asymmetric(curr->index(), (*source)->index())) {
        // a 4; S->L
        std::unique_ptr<Expansion> contrib =
            (*source)->expansion()->S_to_L(curr->center(), (*source)->first(),
                                          (*source)->last());
        curr->expansion()->add_expansion(contrib.get());
      } else {
        // a 1; add it to newcons to push to child of curr, or to handle in the
        // second half of this function.
        newcons.push_back(*source);
      }
    } else {
      if (well_sep_test((*source)->index(), curr->index())) {
        // a 2; M->L
        std::unique_ptr<Expansion> contrib = (*source)->expansion()->M_to_L(
              (*source)->index(), (*source)->size(), curr->index());
        curr->expansion()->add_expansion(contrib.get());
      } else if ((*source)->is_leaf()) {
        newcons.push_back(*source);
      } else {
        // add children for later consideration
        for (size_t i = 0; i < 8; ++i) {
          SourceNode *child = (*source)->child(i);
          if (child != nullptr) {
            newcons.push_back(child);
          }
        }
      }
    }
  }

  if (curr_is_leaf) {
    //If a leaf, we need to also deal with lists 1 and 3
    do {
      std::vector<SourceNode *> updated{};
      for (auto source = newcons.begin(); source != newcons.end(); ++source) {
        if ((*source)->level() < curr->level()) {
          //This would be a larger leaf inherited from the first part of this
          // function.
          assert((*source)->is_leaf());
          (*source)->expansion()->S_to_T((*source)->first(), (*source)->last(),
                                         curr->first(), curr->last());
        } else if (well_sep_test_asymmetric((*source)->index(),
                                            curr->index())) {
          //Otherwise, if the source is well separated, can use M->T; this is
          // list 3
          (*source)->expansion()->M_to_T(curr->first(), curr->last());
        } else if ((*source)->is_leaf()) {
          //Otherwise, if this is a leaf, we must use the direct force
          (*source)->expansion()->S_to_T((*source)->first(), (*source)->last(),
                                         curr->first(), curr->last());
        } else {
          //descend one level and see what is up with the children
          for (size_t child = 0; child < 8; ++child) {
            SourceNode *check = (*source)->child(child);
            if (check != nullptr) {
              updated.push_back(check);
            }
          }
        }
      }
      newcons = std::move(updated);
    } while (newcons.size() != 0);

    curr->expansion()->L_to_T(curr->first(), curr->last());
  } else {
    //Not a leaf, so move newcons into the consider list
    consider = std::move(newcons);
  }
}


bool FMM_Method::refine_test(bool same_sources_and_targets,
    const TargetNode *curr, const std::vector<SourceNode *> &consider) const {
  //
  if (same_sources_and_targets) {
    return true;
  }

  for (auto i = consider.begin(); i != consider.end(); ++i) {
    if ((*i)->level() == curr->level()) {
      if (!well_sep_test((*i)->index(), curr->index()) && !(*i)->is_leaf()) {
        return true;
      }
    }
  }

  return false;
}


bool FMM_Method::well_sep_test_asymmetric(Index smaller, Index larger) const {
  assert(smaller.level() >= larger.level());

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


bool FMM_Method::well_sep_test(Index source, Index target) const {
  assert(source.level() == target.level());

  //When the nodes are the same level, we just need to have at least one
  // index that is different by +/- 2 or more.
  if (abs(source.x() - target.x()) > 1) return true;
  if (abs(source.y() - target.y()) > 1) return true;
  if (abs(source.z() - target.z()) > 1) return true;
  return false;
}


} //namespace dashmm

/// \file src/direct_method.cc
/// \brief Implementation of the direct summation method


#include "include/direct_method.h"


namespace dashmm {


void Direct::generate(SourceNode &curr, const ExpansionRef expand) const {
  //We must at least create an expansion here, otherwise the S->T will not have
  // anything to work from.
  int n_digits = expand.accuracy(); 
  curr.set_expansion(expand.get_new_expansion(Point{0.0, 0.0, 0.0}, n_digits));
}


void Direct::aggregate(SourceNode &curr,
                       const ExpansionRef expand) const {
  //
}

void Direct::inherit(TargetNode &curr, const ExpansionRef expand,
                     size_t which_child) const { 
  int n_digits = expand.accuracy(); 
  curr.set_expansion(expand.get_new_expansion(Point{0.0, 0.0, 0.0}, n_digits));
}

void Direct::process(TargetNode &curr, std::vector<SourceNode> &consider,
                     bool curr_is_leaf) const {
  std::vector<SourceNode> newcons{};
  do {
    for (auto i = consider.begin(); i != consider.end(); ++i) {
      if (i->is_leaf()) {
        if (curr_is_leaf) {
          ExpansionRef expand = i->expansion();
          TargetRef targets = curr.parts();
          SourceRef sources = i->parts();
          expand.S_to_T(sources, targets);
        } else {
          newcons.push_back(*i);
        }
      } else {
        for (size_t j = 0; j < 8; ++j) {
          SourceNode kid = i->child(j);
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


} // namespace dashmm

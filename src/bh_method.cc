/// \file src/bh_method.cc
/// \brief Implementation of BHMethod


#include "include/bh_method.h"


namespace dashmm {


void BH::generate(SourceNode &curr, const ExpansionRef expand) const {
  curr.set_expansion(expand.get_new_expansion(Point{0.0, 0.0, 0.0}));
  ExpansionRef currexp = curr.expansion();
  SourceRef sources = curr.parts();
  double unused = 0.0; 
  currexp.S_to_M(Point{0.0, 0.0, 0.0}, sources, unused);
}


void BH::aggregate(SourceNode &curr, const ExpansionRef expand) const {
  curr.set_expansion(expand.get_new_expansion(Point{0.0, 0.0, 0.0}));
  ExpansionRef currexp = curr.expansion();
  for (size_t i = 0; i < 8; ++i) {
    SourceNode kid = curr.child(i);
    if (kid.is_valid()) {
      ExpansionRef kexp = kid.expansion();  //this data is available in the node
      currexp.M_to_M(kexp, i, 0.0);
    }
  }
}


//NOTE: We pass in curr_is_leaf separate here, as the TargetNode will not have
// been refined yet. Process happens during partition.
void BH::process(TargetNode &curr, std::vector<SourceNode> &consider,
                 bool curr_is_leaf) const {
  std::vector<SourceNode> newcons{};
  double unused = 0.0; 

  do {
    for (auto i = consider.begin(); i != consider.end(); ++i) {
      Point comp_point = nearest(i->center(), curr.center(), curr.size());
      bool can_use = MAC(i->center(), i->size(), comp_point);

      if (can_use) {
        if (!curr_is_leaf) {
          newcons.push_back(*i);
        } else {
          ExpansionRef expand = i->expansion();
          TargetRef targets = curr.parts();
          expand.M_to_T(targets, unused);
        }
      } else if (i->is_leaf()) {
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
          //NOTE: This line will create the SourceNode, and then pull in a
          // local version.
          SourceNode kid = i->child(j);
          if (kid.is_valid()) {
            //NOTE: This will do the same thing. In distrib, that pull might
            // be costly.
            newcons.push_back(kid);
          }
        }
      }
    }

    consider = std::move(newcons);
    newcons.clear();
  } while (curr_is_leaf && !consider.empty());
}


bool BH::MAC(Point exp_point, double size, Point pos) const {
  Point disp = point_sub(pos, exp_point);
  double theta = size / disp.norm();
  return theta < local_->data[0];
}


Point BH::nearest(Point scenter, Point tcenter, double tsize) const {
  Point offset{point_sub(scenter, tcenter)};

  double len{0.5 * tsize};
  double xval{fabs(offset.x()) > len ? len : offset.x()};
  xval += tcenter.x();
  double yval{fabs(offset.y()) > len ? len : offset.y()};
  yval += tcenter.y();
  double zval{fabs(offset.z()) > len ? len : offset.z()};
  zval += tcenter.z();
  return Point{xval, yval, zval};
}


} // namespace dashmm

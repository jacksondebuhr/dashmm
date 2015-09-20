#include <cmath>
#include <vector>
#include <utility>

#include "BH_method.h"
#include "expansion.h"
#include "point.h"


namespace dashmm {


BH_Method::BH_Method(double crit) : crit_{crit} { }


BH_Method::~BH_Method() { }


void BH_Method::generate(SourceNode *curr, const Expansion *expand) const {
  curr->set_expansion(expand->S_to_M(Point{0.0, 0.0, 0.0}, curr->first(),
                                     curr->last()));
}


void BH_Method::aggregate(SourceNode *curr, const Expansion *expand) const {
  curr->set_expansion(expand->get_new_expansion(Point{0.0, 0.0, 0.0}));

  std::vector<std::unique_ptr<Expansion>> exps{};
  std::vector<const Expansion *> refs{};
  for (size_t i = 0; i < 8; ++i) {
    SourceNode *kid = curr->child(i);
    if (kid) {
      exps.push_back(kid->expansion()->M_to_M(i, 0.0));
      refs.push_back(exps[exps.size() - 1].get());
    }
  }
  curr->expansion()->from_sum(refs);
}


void BH_Method::process(TargetNode *curr, std::vector<SourceNode *> &consider,
                        bool curr_is_leaf) const {
  std::vector<SourceNode *> newcons{};

  for (auto i = consider.begin(); i != consider.end(); ++i) {
    Point comp_point = nearest(*i, curr);
    bool can_use = MAC(*i, comp_point);

    if (can_use) {
      //Node can be used
      (*i)->expansion()->M_to_T(curr->first(), curr->last());
    } else if ((*i)->is_leaf()) {
      //Node cannot, and is a leaf
      if (curr_is_leaf) {
        //if curr is a leaf, we should do direct
        (*i)->expansion()->S_to_T((*i)->first(), (*i)->last(),
                                  curr->first(), curr->last());
      } else {
        //otherwise, we put this node in the consider list for curr's children
        // to examine
        newcons.push_back(*i);
      }
    } else {
      //Node cannot, but is not a leaf; so consider its children
      for (size_t j = 0; j < 8; ++j) {
        if ((*i)->child(j)) {
          newcons.push_back((*i)->child(j));
        }
      }
    }
  }

  //then move newcons into consider
  consider = std::move(newcons);

  //If we are at a leaf, the descent needs to finish
  if (curr_is_leaf) {
    finalize(curr, consider);
  }
}


bool BH_Method::MAC(SourceNode *source, Point pos) const {
  Point disp = pos - source->expansion()->center();
  double theta = source->size() / disp.norm();

  return theta < crit_;
}


Point BH_Method::nearest(const SourceNode *source,
                         const TargetNode *target) const {
  Point geocenter{target->center()};
  Point offset{source->expansion()->center() - geocenter};

  //NOTE: This assumes cubic node geometries
  double len{0.5 * (target->high().x() - target->low().x())};
  double xval{fabs(offset.x()) > len ? len : offset.x()};
  xval += geocenter.x();
  double yval{fabs(offset.y()) > len ? len : offset.y()};
  yval += geocenter.y();
  double zval{fabs(offset.z()) > len ? len : offset.z()};
  zval += geocenter.z();
  return Point{xval, yval, zval};
}


void BH_Method::finalize(TargetNode *curr,
                         std::vector<SourceNode *> &consider) const {
  std::vector<SourceNode *> newcons{};

  do {
    for (auto i = consider.begin(); i != consider.end(); ++i) {
      Point comp_point = nearest(*i, curr);
      bool can_use = MAC(*i, comp_point);

      if (can_use) {
        //Node can be used
        (*i)->expansion()->M_to_T(curr->first(), curr->last());
      } else if ((*i)->is_leaf()) {
        //Node cannot, and is a leaf; also, curr must be a leaf as finalize
        // is only called on leaves. So we can do direct here
        (*i)->expansion()->S_to_T((*i)->first(), (*i)->last(),
                                  curr->first(), curr->last());
      } else {
        //Node cannot, but is not a leaf; so consider its children
        for (size_t j = 0; j < 8; ++j) {
          if ((*i)->child(j)) {
            newcons.push_back((*i)->child(j));
          }
        }
      }
    }

    //then move newcons into consider
    consider = std::move(newcons);
    newcons.clear();
  } while (!consider.empty());
}


} //namespace dashmm

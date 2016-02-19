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


#ifndef __DASHMM_BH_METHOD_H__
#define __DASHMM_BH_METHOD_H__


/// \file include/bh_method.h
/// \brief Declaration of BHMethod


#include "include/expansionlco.h"
#include "include/method.h"
#include "include/node.h"


namespace dashmm {


/// A Method to implement classic Barnes-Hut
///
/// It uses the simple critical angle criterion to decide if a given expansion
/// is usable. There is little that needs explanation in this class, as the
/// methods are essentially just those of the abstract base class Method.
template <typename Source, typename Target,
          template <typename, typename> class Expansion>
class BH {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = BH<Source, Target, Expansion>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, BH>;
  using sourcenode_t = SourceNode<Source, Target, Expansion, BH>;
  using targetnode_t = TargetNode<Source, Target, Expansion, BH>;

  BH(double theta) : theta_{theta} {  }

  double theta() const {return theta_;}

  void generate(sourcenode_t &curr, int n_digits) const {
    double scale = 0.0;
    curr.set_expansion(std::unique_ptr<expansion_t>{
        new expansion_t{Point{0.0, 0.0, 0.0}, n_digits}
      });
    expansionlco_t currexp = curr.expansion();
    sourceref_t sources = curr.parts();
    currexp.S_to_M(Point{0.0, 0.0, 0.0}, sources, scale);
  }

  void aggregate(sourcenode_t &curr, int n_digits) const {
    curr.set_expansion(std::unique_ptr<expansion_t>{
        new expansion_t{Point{0.0, 0.0, 0.0}, n_digits}
      });
    expansionlco_t currexp = curr.expansion();
    for (size_t i = 0; i < 8; ++i) {
      sourcenode_t kid = curr.child(i);
      if (kid.is_valid()) {
        expansionlco_t kexp = kid.expansion();
        currexp.M_to_M(kexp, i, 0.0);
      }
    }
  }

  void inherit(targetnode_t &curr, int n_digits, size_t which_child) const {
    curr.set_expansion(std::unique_ptr<expansion_t>{
        new expansion_t{Point{0.0, 0.0, 0.0}, n_digits}
      });
  }

  void process(targetnode_t &curr, std::vector<sourcenode_t> &consider,
               bool curr_is_leaf) const {
    std::vector<sourcenode_t> newcons{ };
    double unused = 0.0;

    do {
      for (auto i = consider.begin(); i != consider.end(); ++i) {
        Point comp_point = nearest(i->center(), curr.center(), curr.size());
        bool can_use = MAC(i->center(), i->size(), comp_point);

        if (can_use) {
          if (!curr_is_leaf) {
            newcons.push_back(*i);
          } else {
            expansionlco_t expand = i->expansion();
            targetlco_t targets = curr.parts();
            expand.M_to_T(targets, unused);
          }
        } else if (i->is_leaf()) {
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
            // NOTE: This line will create the SourceNode, and then pull in a
            // local version.
            sourcenode_t kid = i->child(j);
            if (kid.is_valid()) {
              // NOTE: This will do the same thing. In distrib, that pull might
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

  // BH always calls for refinement
  bool refine_test(bool same_sources_and_targets, const targetnode_t &curr,
                   const std::vector<sourcenode_t> &consider) const {
    return true;
  }

  /// Decide on the usability of an expansion
  ///
  /// This performs the traditional critical angle comparison.
  ///
  /// \param exp_point - the position of the expansion center
  /// \param size - the size of the containing node
  /// \param pos - the position of the point to check for usability.
  ///
  /// \returns - true if the expansion is usable; false otherwise
  bool MAC(Point exp_point, double size, Point pos) const {
    Point disp = point_sub(pos, exp_point);
    double theta = size / disp.norm();
    return theta < local_->data[0];
  }

  /// Finds the point nearest the given point in a given node
  ///
  /// This is used to compute the point in a target node that is nearest
  /// the source expansion center inside a given target node. This allows the
  /// MAC to be used once to verify for all possible points inside the target
  /// node.
  ///
  /// \param scenter - the center of the source expansion
  /// \param tcenter - the center of the target node
  /// \param tsize - the size of the target node
  ///
  /// \returns - the point in the target node closest to the expansion center
  Point nearest(Point scenter, Point tcenter, double tsize) const {
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

 private:
  double theta_;
};


} // namespace dashmm


#endif // __DASHMM_BH_METHOD_H__

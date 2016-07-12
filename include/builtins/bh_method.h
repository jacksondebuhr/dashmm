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


/// \file include/builtins/bh_method.h
/// \brief Declaration of BH method


#include <cstdlib>

#include <vector>

#include "dashmm/arrayref.h"
#include "dashmm/expansionlco.h"
#include "dashmm/point.h"
#include "dashmm/targetlco.h"
#include "dashmm/tree.h"


namespace dashmm {


/// A Method to implement classic Barnes-Hut
///
/// It uses the simple critical angle criterion to decide if a given expansion
/// is usable. If the point of interest is a distance D from the multipole
/// moment center, and the size of the node for which that multipole moment
/// applies is L, then the moments are usable if L / D < theta_c, where
/// theta_c is the critical angle supplied when an instance of this method
/// is constructed.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          typename DistroPolicy = SingleLocality>
class BH {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = BH<Source, Target, Expansion, DistroPolicy>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, BH,
                                      DistroPolicy>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, BH, DistroPolicy>;
  using sourcenode_t = TreeNode<Source, Target, Source, Expansion, BH,
                                DistroPolicy>;
  using targetnode_t = TreeNode<Source, Target, Target, Expansion, BH,
                                DistroPolicy>;
  using sourceref_t = ArrayRef<Source>;

  BH() : theta_{0.0} { }

  /// The BH Method requires a critical angle
  BH(double theta) : theta_{theta} {  }

  /// Return the critical angle
  double theta() const {return theta_;}

  /// In generate, BH will call S->M on the sources in a leaf node.
  void generate(sourcenode_t *curr, DomainGeometry *domain) const {
    curr->dag.add_parts(hpx_get_my_rank()); // TODO: Fix this call to HPX
    curr->dag.add_normal();
    curr->dag.StoM(&curr->dag);
  }

  /// In aggregate, BH will call M->M to combine moments from the children
  /// of the current node.
  void aggregate(sourcenode_t *curr, DomainGeometry *domain) const {
    curr->dag.add_normal();
    for (size_t i = 0; i < 8; ++i) {
      sourcenode_t *kid = curr->child[i];
      if (kid != nullptr) {
        curr->dag.MtoM(&kid->dag);
      }
    }
  }

  /// In inherit, BH does nothing.
  void inherit(targetnode_t *curr, DomainGeometry *domain,
               bool curr_is_leaf) const {
    if (curr_is_leaf) {
      curr->dag.add_parts(hpx_get_my_rank()); // TODO fix this!
    }
  }

  /// In process, BH tests and uses multipole expansions.
  void process(targetnode_t *curr, std::vector<sourcenode_t *> &consider,
               bool curr_is_leaf, DomainGeometry *domain) const {
    std::vector<sourcenode_t *> newcons{ };
    Point ccenter = domain->center_from_index(curr->idx);
    double csize = domain->size_from_level(curr->idx.level());

    do {
      for (auto i = consider.begin(); i != consider.end(); ++i) {
        Point icenter = domain->center_from_index((*i)->idx);
        double isize = domain->size_from_level((*i)->idx.level());
        Point comp_point = nearest(icenter, ccenter, csize);
        bool can_use = MAC(icenter, isize, comp_point);

        if (can_use) {
          if (!curr_is_leaf) {
            newcons.push_back(*i);
          } else {
            (*i)->dag.MtoT(&curr->dag);
          }
        } else if ((*i)->is_leaf()) {
          if (curr_is_leaf) {
            curr->dag.StoT(&(*i)->dag);
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

  // BH always calls for refinement
  bool refine_test(bool same_sources_and_targets, const targetnode_t *curr,
                   const std::vector<sourcenode_t *> &consider) const {
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
    return theta < theta_;
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
  /// The critical angle for this instance of the BH method.
  double theta_;
};


} // namespace dashmm


#endif // __DASHMM_BH_METHOD_H__

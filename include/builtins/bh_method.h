// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
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
#include "dashmm/defaultpolicy.h"
#include "dashmm/expansionlco.h"
#include "dashmm/point.h"
#include "dashmm/targetlco.h"
#include "dashmm/tree.h"

#include "builtins/bhdistro.h"


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
          template <typename, typename> class Expansion>
class BH {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = BH<Source, Target, Expansion>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, BH>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, BH>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using sourceref_t = ArrayRef<Source>;

  using distropolicy_t = BHDistro;

  BH() : theta_{0.0} { }

  /// The BH Method requires a critical angle
  BH(double theta) : theta_{theta} {  }

  /// Return the critical angle
  double theta() const {return theta_;}

  /// In generate, BH will call S->M on the sources in a leaf node.
  void generate(sourcenode_t *curr, DomainGeometry *domain) const {
    curr->dag.add_parts();
    assert(curr->dag.add_normal() == true);
    curr->dag.StoM(&curr->dag,
                   expansion_t::weight_estimate(Operation::StoM));
  }

  /// In aggregate, BH will call M->M to combine moments from the children
  /// of the current node.
  void aggregate(sourcenode_t *curr, DomainGeometry *domain) const {
    assert(curr->dag.add_normal() == true);
    for (size_t i = 0; i < 8; ++i) {
      sourcenode_t *kid = curr->child[i];
      if (kid != nullptr) {
        curr->dag.MtoM(&kid->dag,
                       expansion_t::weight_estimate(Operation::MtoM));
      }
    }
  }

  /// In inherit, BH does nothing.
  void inherit(targetnode_t *curr, DomainGeometry *domain,
               bool curr_is_leaf) const {
    if (curr_is_leaf) {
      curr->dag.add_parts();
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
            (*i)->dag.MtoT(&curr->dag,
                           expansion_t::weight_estimate(Operation::MtoT));
          }
        } else if ((*i)->is_leaf()) {
          if (curr_is_leaf) {
            curr->dag.StoT(&(*i)->dag,
                           expansion_t::weight_estimate(Operation::StoT));
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

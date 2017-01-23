// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_DOMAIN_GEOMETRY_H__
#define __DASHMM_DOMAIN_GEOMETRY_H__


/// \file
/// \brief Object to represent the top-level domain


#include "dashmm/index.h"
#include "dashmm/point.h"


namespace dashmm {


/// The geometry of the overall domain covered by the source and target points
///
/// This is used to specify the geometry of the root node for both the souce
/// and target trees.
class DomainGeometry {
 public:
  /// Default constructor
  DomainGeometry() : low_{Point{0.0, 0.0, 0.0}}, size_{0.0} { }

  /// Constructor from low corner and side length
  DomainGeometry(Point low, double size)
      : low_{low}, size_{size} { }

  /// Cubifying constructor
  ///
  /// Given the low and high corners of a rectangular region, this constructor
  /// will initialize the DomainGeometry to be the cube that shares a center
  /// with the specified rectangular region, but which contains the entire
  /// rectangle. Optionally, one can increase the size of the resulting domain
  /// by some factor.
  ///
  /// \param low_rect - the low corner of the rectangular region
  /// \param high_rect - the high corner of the rectangular region
  /// \param factor - factor by which to enlarge the resulting cube
  DomainGeometry(Point low_rect, Point high_rect, double factor = 1.0);

  /// The side length of the cube represented by this object
  double size() const {return size_;}

  /// The low corner of the cube represented by this object
  Point low() const {return low_;}

  /// The high corner of the cube represented by this object
  Point high() const {
    return Point{low_.x() + size_, low_.y() + size_, low_.z() + size_};
  }

  /// The center of the cube represented by this object
  Point center() const {
    return Point{low_.x() + 0.5 * size_, low_.y() + 0.5 * size_,
                 low_.z() + 0.5 * size_};
  }

  /// Return the low point of a subdivision of this domain.
  ///
  /// This computes the low corner of a given partition of this domain,
  /// specified as an Index.
  ///
  /// \param idx - the index of the subdivision.
  ///
  /// \returns - the location of the low corner of the specified subdivision.
  Point low_from_index(Index idx) const;

  /// Return the high point of a subdivision of this domain.
  ///
  /// This computes the high corner of a given partition of this domain,
  /// specified as an Index.
  ///
  /// \param idx - the index of the subdivision.
  ///
  /// \returns - the location of the high corner of the specified subdivision.
  Point high_from_index(Index idx) const;

  /// Return the center point of a subdivision of this domain.
  ///
  /// This computes the center of a given partition of this domain,
  /// specified as an Index.
  ///
  /// \param idx - the index of the subdivision.
  ///
  /// \returns - the location of the center of the specified subdivision.
  Point center_from_index(Index idx) const;

  /// Return the size of a subdivision given the level
  ///
  /// \param level - the level in question
  ///
  /// \returns - the size of nodes at that level.
  double size_from_level(int level) const;

 private:
  Point low_;
  double size_;
};


} // namespace dashmm


#endif

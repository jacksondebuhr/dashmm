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


/// \file src/domaingeometry.cc
/// \brief Implementation of DomainGeometry


#include <cmath>

#include "dashmm/domaingeometry.h"


namespace dashmm {


DomainGeometry::DomainGeometry(Point low_rect, Point high_rect, double factor)
    : low_{0.0, 0.0, 0.0}, size_{0.0} {
  Point center = point_add(low_rect.scale(0.5), high_rect.scale(0.5));
  Point sizes = point_sub(high_rect, low_rect);
  double max_size = sizes.x() > sizes.y() ? sizes.x() : sizes.y();
  max_size = max_size > sizes.z() ? max_size : sizes.z();
  max_size *= factor;
  size_ = max_size;
  max_size *= 0.5;
  low_ = Point{center.x() - max_size, center.y() - max_size,
               center.z() - max_size};
}


Point DomainGeometry::low_from_index(Index idx) const {
  double level_size = size_from_level(idx.level());
  return Point{idx.x() * level_size + low_.x(),
               idx.y() * level_size + low_.y(),
               idx.z() * level_size + low_.z()};
}


Point DomainGeometry::high_from_index(Index idx) const {
  double level_size = size_from_level(idx.level());
  return Point{(idx.x() + 1) * level_size + low_.x(),
               (idx.y() + 1) * level_size + low_.y(),
               (idx.z() + 1) * level_size + low_.z()};
}


Point DomainGeometry::center_from_index(Index idx) const {
  double level_size = size_from_level(idx.level());
  return Point{(idx.x() + 0.5) * level_size + low_.x(),
               (idx.y() + 0.5) * level_size + low_.y(),
               (idx.z() + 0.5) * level_size + low_.z()};
}


double DomainGeometry::size_from_level(int level) const {
  return size_ / static_cast<double>(1 << level);
}


} // namespace dashmm

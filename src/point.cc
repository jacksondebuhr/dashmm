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


/// \file src/point.cc
/// \brief Implemention of Point operations


#include "dashmm/point.h"


namespace dashmm {


double point_dot(const Point &left, const Point &right) {
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}


Point point_add(const Point &left, const Point &right) {
  return Point{left[0] + right[0], left[1] + right[1], left[2] + right[2]};
}


Point point_sub(const Point &left, const Point &right) {
  return Point{left[0] - right[0], left[1] - right[1], left[2] - right[2]};
}


} // namespace dashmm

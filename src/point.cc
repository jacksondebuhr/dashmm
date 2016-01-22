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


/// \file src/point.cc
/// \brief Implemention of Point operations


#include "include/point.h"


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

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

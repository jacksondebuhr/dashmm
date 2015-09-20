#include "point.h"


namespace dashmm {

double operator*(const Point &left, const Point &right) {
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

Point operator+(const Point &left, const Point &right) {
  return Point{left[0] + right[0], left[1] + right[1], left[2] + right[2]};
}

Point operator-(const Point &left, const Point &right) {
  return Point{left[0] - right[0], left[1] - right[1], left[2] - right[2]};
}


} // namespace dashmm

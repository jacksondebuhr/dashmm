#ifndef __DASHMM_POINT_H__
#define __DASHMM_POINT_H__


#include <cmath>


namespace dashmm {


class Point {
 public:
  Point(double x, double y, double z) : pos_{x, y, z} { }
  Point(const Point &pt) : pos_{pt.x(), pt.y(), pt.z()} { }

  Point scale(double c) const {
    return Point{c * pos_[0], c * pos_[1], c * pos_[2]};
  }
  double operator[](size_t i) const {return pos_[i];}
  double x() const {return pos_[0];}
  double y() const {return pos_[1];}
  double z() const {return pos_[2];}
  double norm() const {
    return sqrt(pos_[0] * pos_[0] + pos_[1] * pos_[1] + pos_[2] * pos_[2]);
  }

 private:
  double pos_[3];
};


double point_dot(const Point &left, const Point &right);

Point point_add(const Point &left, const Point &right);

Point point_sub(const Point &left, const Point &right);


} //namespace dashmm


#endif // __DASHMM_POINT_H__

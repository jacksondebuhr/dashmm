#ifndef __DASHMM_DOMAIN_GEOMETRY_H__
#define __DASHMM_DOMAIN_GEOMETRY_H__


#include "include/point.h"


namespace dashmm {


class DomainGeometry {
 public:
  DomainGeometry(Point low, double size)
      : low_{low}, size_{size} { }
  //This one takes in rectangular bounds and makes a cube from them
  DomainGeometry(Point low_rect, Point high_rect);

  double size() const {return size_;}
  Point low() const {return low_;}
  Point high() const {
    return Point{low_.x() + size_, low_.y() + size_, low_.z() + size_};
  }
  Point center() const {
    return Point{low_.x() + size_, low_.y() + size_, low_.z() + size_};
  }

  Point low_from_index(int ix, int iy, int iz, int level) const;
  Point high_from_index(int ix, int iy, int iz, int level) const;
  Point center_from_index(int ix, int iy, int iz, int level) const;
  double size_from_level(int level) const;

 private:
  Point low_;
  double size_;
};


} // namespace dashmm


#endif

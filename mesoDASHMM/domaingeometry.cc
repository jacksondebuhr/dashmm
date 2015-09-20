#include "domaingeometry.h"

#include <cmath>


namespace dashmm {


DomainGeometry::DomainGeometry(Point low_rect, Point high_rect)
    : low_{0.0, 0.0, 0.0}, size_{0.0} {
  Point center{low_rect.x() * 0.5 + high_rect.x() * 0.5,
               low_rect.y() * 0.5 + high_rect.y() * 0.5,
               low_rect.z() * 0.5 + high_rect.z() * 0.5};
  Point diff{high_rect.x() - low_rect.x(), high_rect.y() - low_rect.y(),
             high_rect.z() - low_rect.z()};
  double max_size = diff.x() > diff.y() ? diff.x() : diff.y();
  max_size = max_size > diff.z() ? max_size : diff.z();
  size_ = max_size;
  max_size *= 0.5;
  low_ = Point{center.x() - max_size, center.y() - max_size,
               center.z() - max_size};
}


Point DomainGeometry::low_from_index(int ix, int iy, int iz, int level) const {
  double level_size = size_from_level(level);
  return Point{ix * level_size + low_.x(), 
      iy * level_size + low_.y(), 
      iz * level_size + low_.z()};
}


Point DomainGeometry::high_from_index(int ix, int iy, int iz, int level) const {
  double level_size = size_from_level(level);
  return Point{(ix + 1) * level_size + low_.x(), 
      (iy + 1) * level_size + low_.y(),
      (iz + 1) * level_size + low_.z()};
}


Point DomainGeometry::center_from_index(int ix, int iy, int iz,
                                        int level) const {
  double level_size = size_from_level(level);
  return Point{(static_cast<double>(ix) + 0.5) * level_size + low_.x(),
               (static_cast<double>(iy) + 0.5) * level_size + low_.y(),
               (static_cast<double>(iz) + 0.5) * level_size + low_.z()};
}


double DomainGeometry::size_from_level(int level) const {
  return size_ / static_cast<double>(1 << level);
}


} //namespace dashmm

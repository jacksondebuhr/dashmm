#include "include/domaingeometry.h"

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

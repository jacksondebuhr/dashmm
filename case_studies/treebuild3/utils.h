#ifndef __TREBUILD_3_UTILS_H__
#define __TREBUILD_3_UTILS_H__


#include <cmath>
#include <cstdlib>
#include <cstdint>

#include "dashmm/point.h"
using dashmm::Point;

#include "dashmm/index.h"
using dashmm::Index;

#include "dashmm/array.h"


void set_point_in_cube(Point &p);
void set_point_on_sphere(Point &p);

uint64_t split(unsigned k);
uint64_t morton_key(unsigned x, unsigned y, unsigned z);

dashmm::Array<Point> generate_points(char scaling, char datatype, int nsrc,
                                      int nseed, int shift);


#endif // __TREEBUILD_3_UTILS_H__
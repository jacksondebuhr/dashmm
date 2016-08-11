#ifndef __TREBUILD_3_UTILS_H__
#define __TREBUILD_3_UTILS_H__


#include "dashmm/point.h"
#include "dashmm/array.h"


void set_point_in_cube(dashmm::Point &p);
void set_point_on_sphere(dashmm::Point &p);
dashmm::Array<dashmm::Point> generate_points(char scaling, char datatype,
                                             int nsrc, int nseed, int shift);


#endif // __TREEBUILD_3_UTILS_H__

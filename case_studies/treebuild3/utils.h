#ifndef __TREBUILD_3_UTILS_H__
#define __TREBUILD_3_UTILS_H__


#include <cmath>
#include <cstdlib>
#include <cstdint>

#include "dashmm/point.h"
using dashmm::Point;

#include "dashmm/index.h"
using dashmm::Index;


/// Meta data for the array.
///
/// The array is a cyclic allocation of ArrayMetaData objects, one per locality
/// that inside point to the local data (data). This gives the size and the
/// count of the records as well as the pointer to local memory.
struct ArrayMetaData {
  size_t size;
  int count;
  char *data;
};


void set_point_in_cube(Point &p);
void set_point_on_sphere(Point &p);

uint64_t split(unsigned k);
uint64_t morton_key(unsigned x, unsigned y, unsigned z);

void generate_weak_scaling_input(ArrayMetaData *meta_s, int nsrc_per_rank,
                                 ArrayMetaData *meta_t, int ntar_per_rank,
                                 char datatype, int seed);
void generate_strong_scaling_input(ArrayMetaData *meta_s, int nsrc,
                                   ArrayMetaData *meta_t, int ntar,
                                   char datatype, int rank, int num_ranks,
                                   int nseed);


#endif // __TREEBUILD_3_UTILS_H__
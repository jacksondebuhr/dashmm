#include "utils.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>

using dashmm::Point;
using dashmm::Array;
using dashmm::ArrayRef;
using dashmm::ArrayData;


void set_point_in_cube(Point &p) {
  double x = 1.0 * rand() / RAND_MAX - 0.5;
  double y = 1.0 * rand() / RAND_MAX - 0.5;
  double z = 1.0 * rand() / RAND_MAX - 0.5;
  p = Point{x, y, z};
}

void set_point_on_sphere(Point &p) {
  double theta = 1.0 * rand() / RAND_MAX * M_PI_2;
  double phi = 1.0 * rand() / RAND_MAX * M_PI * 2;
  double x = sin(theta) * cos(phi);
  double y = sin(theta) * sin(phi);
  double z = cos(theta);
  p = Point{x, y, z};
}

Point *generate_weak_scaling_input(int n, char datatype, int seed) {
  Point *retval = new Point[n]();

  srand(seed);
  std::function<void(Point &)> set_point =
    (datatype == 'c' ? set_point_in_cube : set_point_on_sphere);

  for (int i = 0; i < n; ++i) {
    set_point(retval[i]);
  }

  return retval;
}

Point *generate_strong_scaling_input(int n, char datatype, int rank,
                                     int num_ranks, int nseed, int shift) {
  // Divide nsrc and ntar into nseed chunks. Use seed (2 + r) to populate
  // the points in each chunk, where 0 <= r < nseed

  std::function<void(Point &)> set_point =
    (datatype == 'c' ? set_point_in_cube : set_point_on_sphere);

  int q1, r1, q2, r2, nsrc_curr_rank, s1, s2;

  // Generate source points
  // nsrc = q1 * num_ranks + r1 = q2 * nseed + r2
  q1 = n / num_ranks;
  r1 = n % num_ranks;
  q2 = n / nseed;
  r2 = n % nseed;
  nsrc_curr_rank = (rank < r1 ? q1 + 1 : q1);
  s1 = (rank < r1 ? (q1 + 1) * rank : (q1 + 1) * r1 + q1 * (rank - r1));
  s2 = s1 + nsrc_curr_rank - 1;

  Point *retval = new Point[nsrc_curr_rank]();

  for (int r = 0; r < nseed; ++r) {
    int nsrc_curr_seed = (r < r2 ? q2 + 1 : q2);
    int s3 = (r < r2 ? (q2 + 1) * r : (q2 + 1) * r2 + q2 * (r - r2));
    int s4 = s3 + nsrc_curr_seed - 1;
    if (s4 < s1) {
      continue;
    } else if (s3 > s2) {
      break;
    } else {
      srand(r + 2 + shift);
      for (int s = s3; s <= s4; ++s) {
        if (s >= s1 && s <= s2) {
          set_point(retval[s - s1]);
        }
      }
    }
  }

  return retval;
}

int point_count(char scaling, int n) {
  if (scaling == 'w') {
    return n;
  } else {
    int my_rank = hpx_get_my_rank();
    int num_ranks = hpx_get_num_ranks();
    int nper = n / num_ranks;
    int remain = n % num_ranks;
    return (my_rank < remain ? nper + 1 : nper);
  }
}

Array<Point> generate_points(char scaling, char datatype, int nsrc,
                                      int nseed, int shift) {
  size_t my_nsrc = point_count(scaling, nsrc);

  Array<Point> retval{};
  retval.allocate(my_nsrc);

  Point *data{nullptr};
  if (scaling == 'w') {
    data = generate_weak_scaling_input(my_nsrc, datatype,
                                       hpx_get_my_rank() + 2);
  } else {
    data = generate_strong_scaling_input(nsrc, datatype, hpx_get_my_rank(),
                                         hpx_get_num_ranks(), nseed, shift);
  }

  retval.put(0, my_nsrc, data);

  delete [] data;

  return retval;
}

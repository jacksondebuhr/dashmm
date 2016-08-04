#include "utils.h"

#include <functional>

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

uint64_t split(unsigned k) {
  uint64_t split = k & 0x1fffff;
  split = (split | split << 32) & 0x1f00000000ffff;
  split = (split | split << 16) & 0x1f0000ff0000ff;
  split = (split | split << 8)  & 0x100f00f00f00f00f;
  split = (split | split << 4)  & 0x10c30c30c30c30c3;
  split = (split | split << 2)  & 0x1249249249249249;
  return split;
}

uint64_t morton_key(unsigned x, unsigned y, unsigned z) {
  uint64_t key = 0;
  key |= split(x) | split(y) << 1 | split(z) << 2;
  return key;
}

void generate_weak_scaling_input(ArrayMetaData *meta_s, int nsrc_per_rank,
                                 ArrayMetaData *meta_t, int ntar_per_rank,
                                 char datatype, int seed) {
  meta_s->size = sizeof(Point);
  meta_s->count = nsrc_per_rank;
  meta_s->data = new char[sizeof(Point) * nsrc_per_rank]();

  meta_t->size = sizeof(Point);
  meta_t->count = ntar_per_rank;
  meta_t->data = new char[sizeof(Point) * ntar_per_rank]();

  Point *p_s = reinterpret_cast<Point *>(meta_s->data);
  Point *p_t = reinterpret_cast<Point *>(meta_t->data);

  srand(seed);
  std::function<void(Point &)> set_point =
    (datatype == 'c' ? set_point_in_cube : set_point_on_sphere);

  for (int i = 0; i < nsrc_per_rank; ++i) {
    set_point(p_s[i]);
  }

  for (int i = 0; i < ntar_per_rank; ++i) {
    set_point(p_t[i]);
  }
}

void generate_strong_scaling_input(ArrayMetaData *meta_s, int nsrc,
                                   ArrayMetaData *meta_t, int ntar,
                                   char datatype, int rank, int num_ranks,
                                   int nseed) {
  // Divide nsrc and ntar into nseed chunks. Use seed (2 + r) to populate
  // the points in each chunk, where 0 <= r < nseed

  std::function<void(Point &)> set_point =
    (datatype == 'c' ? set_point_in_cube : set_point_on_sphere);

  int q1, r1, q2, r2, nsrc_curr_rank, s1, s2, ntar_curr_rank, t1, t2;

  // Generate source points
  // nsrc = q1 * num_ranks + r1 = q2 * nseed + r2
  q1 = nsrc / num_ranks;
  r1 = nsrc % num_ranks;
  q2 = nsrc / nseed;
  r2 = nsrc % nseed;
  nsrc_curr_rank = (rank < r1 ? q1 + 1 : q1);
  s1 = (rank < r1 ? (q1 + 1) * rank : (q1 + 1) * r1 + q1 * (rank - r1));
  s2 = s1 + nsrc_curr_rank - 1;

  meta_s->size = sizeof(Point);
  meta_s->count = nsrc_curr_rank;
  meta_s->data = new char[sizeof(Point) * nsrc_curr_rank]();
  Point *p_s = reinterpret_cast<Point *>(meta_s->data);

  for (int r = 0; r < nseed; ++r) {
    int nsrc_curr_seed = (r < r2 ? q2 + 1 : q2);
    int s3 = (r < r2 ? (q2 + 1) * r : (q2 + 1) * r2 + q2 * (r - r2));
    int s4 = s3 + nsrc_curr_seed - 1;
    if (s4 < s1) {
      continue;
    } else if (s3 > s2) {
      break;
    } else {
      srand(r + 2);
      for (int s = s3; s <= s4; ++s) {
        if (s >= s1 && s <= s2)
          set_point(p_s[s - s1]);
      }
    }
  }

  // Generate target points
  // ntar = q1 * num_ranks + r1 = q2 * nseed + r2
  q1 = ntar / num_ranks;
  r1 = ntar % num_ranks;
  q2 = ntar / nseed;
  r2 = ntar % nseed;
  ntar_curr_rank = (rank < r1 ? q1 + 1 : q1);
  t1 = (rank < r1 ? (q1 + 1) * rank : (q1 + 1) * r1 + q1 * (rank - r1));
  t2 = t1 + ntar_curr_rank - 1;

  meta_t->size = sizeof(Point);
  meta_t->count = ntar_curr_rank;
  meta_t->data = new char[sizeof(Point) * ntar_curr_rank]();
  Point *p_t = reinterpret_cast<Point *>(meta_t->data);

  for (int r = 0; r < nseed; ++r) {
    int ntar_curr_seed = (r < r2 ? q2 + 1 : q2);
    int t3 = (r < r2 ? (q2 + 1) * r : (q2 + 1) * r2 + q2 * (r - r2));
    int t4 = t3 + ntar_curr_seed - 1;
    if (t4 < t1) {
      continue;
    } else if (t3 > t2) {
      break;
    } else {
      srand(r + 2);
      for (int t = t3; t <= t4; ++t) {
        if (t >= t1 && t <= t2)
          set_point(p_t[t - t1]);
      }
    }
  }
}

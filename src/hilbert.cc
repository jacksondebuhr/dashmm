// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================

#include "dashmm/hilbert.h"

#include <cassert>

#include <algorithm>
#include <array>
#include <vector>
#include <utility>

namespace dashmm {

namespace {

  // Group transformations for hilbert subcubes
  constexpr int h_indicies[8] {0, 1, 3, 2, 7, 6, 4, 5};
  constexpr int h_transforms[8][8] {
    {0, 7, 4, 3, 2, 5, 6, 1},
    {0, 3, 4, 7, 6, 5, 2, 1},
    {0, 3, 4, 7, 6, 5, 2, 1},
    {2, 1, 0, 3, 4, 7, 6, 5},
    {2, 1, 0, 3, 4, 7, 6, 5},
    {6, 5, 2, 1, 0, 3, 4, 7},
    {6, 5, 2, 1, 0, 3, 4, 7},
    {6, 1, 2, 5, 4, 3, 0, 7}
  };

  // Compute the Hilbert key for a given set of indicies
  int hilbert_key(unsigned x, unsigned y, unsigned z, unsigned slen) {
    std::array<int, 8> offsets{0, 1, 2, 3, 4, 5, 6, 7};
    std::array<int, 8> offsets_next{};
    int key = 0;

    for (int s = slen; s > 1; s /= 2) {
      int m_index = 4 * (z % s >= s / 2)
                    + 2 * (y % s >= s / 2)
                    + (x % s >= s / 2);
      key += (s * s * s) / 8 * offsets[h_indicies[m_index]];
      for (int i = 0; i < 8; ++i) {
        offsets_next[i] =
            h_transforms[offsets[h_indicies[m_index]]][offsets[i]];
      }
      std::swap(offsets_next, offsets);
    }

    return key;
  }


  std::vector<int> partition(const std::vector<int> &counts,
                             int np,
                             int len) {
    std::vector<int> retval(len);

    int *cumulative = new int [len + 1];
    cumulative[0] = 0;
    for (int i = 1; i <= len; ++i) {
      cumulative[i] = counts[i - 1] + cumulative[i - 1];
    }

    int target_share = cumulative[len] / np;

    std::vector<int> cuts(np);
    for (int split = 1; split < np; ++split) {
      int my_target = target_share * split;
      auto fcn = [&my_target] (const int &a) -> bool {return a < my_target;};
      int *splitter = std::partition_point(cumulative, &cumulative[len + 1],
                                           fcn);
      int upper_delta = *splitter - my_target;
      assert(upper_delta >= 0);
      int lower_delta = my_target - *(splitter - 1);
      assert(lower_delta >= 0);
      int split_index = upper_delta < lower_delta
                          ? splitter - cumulative
                          : (splitter - cumulative) - 1;
      --split_index;
      assert(split_index >= 0);

      cuts[split - 1] = split_index;
    }
    cuts[np - 1] = len - 1;

    delete [] cumulative;

    // use cut points to set return values
    for (int i = 0; i < np; ++i) {
      int b = i == 0 ? 0 : cuts[i - 1] + 1;
      int e = cuts[i];
      for (int j = b; j <= e; ++j) {
        retval[j] = i;
      }
    }

    return retval;
  }

} // {anonymous}


int *distribute_points_hilbert(int num_ranks,
                               const int *global,
                               int len,
                               int lvl) {
  int klen = 1 << lvl;  // the number of nodes in one direction is 2^lvl

  const int *s = global; // Source counts
  const int *t = &global[len]; // Target counts

  // accumulate the counts and setup an inverse mapping
  std::vector<int> counts(len);
  std::vector<int> origin(len);
  for (int zidx = 0; zidx < klen; ++zidx) {
    int zpart = zidx * klen * klen;
    for (int yidx = 0; yidx < klen; ++yidx) {
      int ypart = zpart + yidx * klen;
      for (int xidx = 0; xidx < klen; ++xidx) {
        int i = xidx + ypart;
        int hidx = hilbert_key(xidx, yidx, zidx, klen);
        counts[hidx] = s[i] + t[i];
        origin[hidx] = i;
      }
    }
  }

  // Break up the segments
  std::vector<int> rm = partition(counts, num_ranks, len);

  // Map back into the original index
  int *retval = new int[len];
  for (int hidx = 0; hidx < len; ++hidx) {
    retval[origin[hidx]] = rm[hidx];
  }

  return retval;
}


} // dashmm

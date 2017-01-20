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


#ifndef __DASHMM_BUILTIN_TABLE_H__
#define __DASHMM_BUILTIN_TABLE_H__


/// \file
/// \brief Declaration of comparator for builtin tables


namespace dashmm {


struct builtin_cmp {
  bool operator()(const double &a, const double &b) const {
    // For builtin kernels, the smallest gap between key values of the
    // rotation matrix map is 0.01. This operator compares the first 6
    // digits and should be sufficient.
    double aa = floor(a * 1000000.0) / 100000.0;
    double bb = floor(b * 1000000.0) / 100000.0;
    if (aa == bb)
      return false;
    return aa < bb;
  }
};

using builtin_map_t = std::map<double, double *, builtin_cmp>;


} // namespace dashmm


#endif // __DASHMM_BUILTIN_TABLE_H__

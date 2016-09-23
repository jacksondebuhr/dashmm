// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_MERGE_SHIFT_H__
#define __DASHMM_MERGE_SHIFT_H__


/// \file
/// \brief Declaration of list enum used in FMM97 method


enum MergedList {
  uall = 0, ///< +z direction list for all boxes
  u1234 = 1, ///< +z direction list for boxes 1, 2, 3, 4
  nall = 2, ///< +y direction list for all boxes
  n1256 = 3, ///< +y direction list for boxes 1, 2, 5, 6
  n12 = 4, ///< +y direction list for boxes 1, 2
  n56 = 5, ///< +y direction list for boxes 5, 6
  eall = 6, ///< +x direction list for all boxes
  e1357 = 7, ///< +x direction list for boxes 1, 3, 5, 7
  e13 = 8, ///< +x direction list for boxes 1, 3
  e57 = 9, ///< +x direction list for boxes 5, 7
  e1 = 10, ///< +x direction list for box 1
  e3 = 11, ///< +x direction list for box 3
  e5 = 12, ///< +x direction list for box 5
  e7 = 13, ///< +x direction list for box 7
  dall = 14, ///< -z direction list for all boxes
  d5678 = 15, ///< -z direction list for boxes 5, 6, 7, 8
  sall = 16, ///< -y direction list for all boxes
  s3478 = 17, ///< -y direction list for boxes 3, 4, 7, 8
  s34 = 18, ///< -y direction list for boxes 3, 4
  s78 = 19, ///< -y direction list for boxes 7, 8
  wall = 20, ///< -x direction list for all boxes
  w2468 = 21, ///< -x direction list for boxes 2, 4, 6, 8
  w24 = 22, ///< -x direction list for boxes 2, 4
  w68 = 23, ///< -x direction list for boxes 6, 8
  w2 = 24, ///< -x direction list for box 2
  w4 = 25, ///< -x direction list for box 4
  w6 = 26, ///< -x direction list for box 6
  w8 = 27, ///< -x direction list for box 8
};

extern const int merge_and_shift_table[6][6][6][3];

#endif // __DASHMM_MERGE_SHIFT_H__

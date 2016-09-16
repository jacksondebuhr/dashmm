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


#ifndef __DASHMM_INDEX_H__
#define __DASHMM_INDEX_H__


/// \file include/dashmm/index.h
/// \brief An index type for specifying relative node geometry


namespace dashmm {


/// An Index for specifying the geometry of a node exactly using integers
///
/// This object specifies the level and the location on that level of a
/// node. Indices require 3 components in three dimensions. At level 0,
/// the only possible index is (0, 0, 0). At level 1, there are eight
/// possible indices, with each component having a value of 0 or 1.
///
/// To convert the index into real positions, a DomainGeometry object is
/// needed. See the DomainGeometry documentation for more.
class Index {
 public:
  /// Construct the index from components and level
  Index(int ix = 0, int iy = 0, int iz = 0, int level = 0)
      : idx_{ix, iy, iz}, level_{level} { }

  /// Return the x index.
  int x() const {return idx_[0];}

  /// Return the y index.
  int y() const {return idx_[1];}

  /// Return the z index.
  int z() const {return idx_[2];}

  /// Return the level of the index.
  int level() const {return level_;}

  /// Compute the Index of this Index's parent
  ///
  /// This will compute the index of one of the ancestors of this index.
  ///
  /// \param num - the number of steps upward in the hierarchy to compute.
  ///
  /// \returns - the index which is the \param num parent of this index.
  Index parent(int num = 1) const {
    return Index{idx_[0] >> num, idx_[1] >> num, idx_[2] >> num,
                 level_ > num ? level_ - num : 0};
  }

  /// Compute the child Index relative to this Index.
  ///
  /// Given which child offset from this Index, this computes the Index of
  /// the child. The ordering of the children of a node are as follows,
  /// where the triple of [+-] indicates which half of the node in each
  /// direction the given child occupies:
  /// 0: ---; 1: +--; 2: -+-; 3: ++-;
  /// 4: --+; 5: +-+; 6: -++; 7: +++;
  ///
  /// \param which - the child Index to compute
  ///
  /// \returns - the Index of the child
  Index child(int which) const {
    return Index{(idx_[0] << 1) + (which & 1 ? 1 : 0),
                 (idx_[1] << 1) + (which & 2 ? 1 : 0),
                 (idx_[2] << 1) + (which & 4 ? 1 : 0),
                 level_ + 1};
  }

  /// Compute which child of this Index's parent this index is
  ///
  /// The ordering of the children of a node are as follows,
  /// where the triple of [+-] indicates which half of the node in each
  /// direction the given child occupies:
  /// 0: ---; 1: +--; 2: -+-; 3: ++-;
  /// 4: --+; 5: +-+; 6: -++; 7: +++;
  ///
  /// \returns - which child this Index is
  int which_child() const {
    int xval = idx_[0] % 2;
    int yval = (idx_[1] % 2) << 1;
    int zval = (idx_[2] % 2) << 2;
    return (xval + yval + zval);
  }

  /// Equality operator
  bool operator==(const Index &other) {
    return (level_ == other.level_
              && idx_[0] == other.idx_[0]
              && idx_[1] == other.idx_[1]
              && idx_[2] == other.idx_[2]);
  }

 private:
  int idx_[3];
  int level_;
};


} // namespace dashmm


#endif // __DASHMM_INDEX_H__

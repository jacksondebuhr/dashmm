// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_INDEX_H__
#define __DASHMM_INDEX_H__


/// \file include/index.h
/// \brief An index type for specifying relative node geometry


namespace dashmm {


/// An Index for specifying the geometry of a node exactly using integers
///
/// This object specifies the level and the location on that level of a
/// node. Indices require 3 components in three dimensions. At level 0,
/// the only possible index is (0, 0, 0). At level 1, there are eight
/// possible indices, with each component having a value of 0 or 1.
class Index {
 public:
  /// construct the index from components and level
  Index(int ix, int iy, int iz, int level)
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

 private:
  int idx_[3];
  int level_;
};


} // namespace dashmm

#endif // __DASHMM_INDEX_H__

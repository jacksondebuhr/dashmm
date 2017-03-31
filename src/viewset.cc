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


#include "dashmm/viewset.h"


/// \file
/// \brief Implementation of ViewSet


namespace dashmm {


void ViewSet::clear() {
  views_.clear();
  role_ = kNoRoleNeeded;
}


void ViewSet::add_view(int index) {
  add_view(index, 0, nullptr);
}


void ViewSet::add_view(int index, size_t bytes, char *data) {
  views_.push_back(View{index, bytes, data});
}


size_t ViewSet::bytes() const {
  size_t retval = 0;
  for (size_t i = 0; i < views_.size(); ++i) {
    retval += views_[i].bytes;
  }

  retval += sizeof(ExpansionRole);          // the role
  retval += sizeof(int);                    // for count
  retval += 4 * sizeof(double);             // for the center and the scale
  retval += views_.size() * sizeof(int);
  retval += views_.size() * sizeof(size_t);

  return retval;
}


void ViewSet::serialize(WriteBuffer &buffer) {
  bool e = buffer.write(role_);
  assert(e);
  e = buffer.write(count());
  assert(e);
  double oval = center_.x();
  e = buffer.write(oval);
  assert(e);
  oval = center_.y();
  e = buffer.write(oval);
  assert(e);
  oval = center_.z();
  e = buffer.write(oval);
  assert(e);
  e = buffer.write(scale_);
  assert(e);

  // then the view index and size for each
  for (int i = 0; i < count(); ++i) {
    e = buffer.write(view_index(i));
    assert(e);
    e = buffer.write(view_bytes(i));
    assert(e);
  }

  // then the data
  for (int i = 0; i < count(); ++i) {
    e = buffer.write(view_data(i), view_bytes(i));
    assert(e);
  }
}


void ViewSet::interpret(ReadBuffer &buffer) {
  // NOTE: this clears out the ViewSet
  views_.clear();

  bool e = buffer.read(&role_);
  assert(e);

  int ct{0};
  e = buffer.read((char *)&ct, sizeof(ct));
  assert(e);

  double ival[3];
  e = buffer.read(&ival[0]);
  assert(e);
  e = buffer.read(&ival[1]);
  assert(e);
  e = buffer.read(&ival[2]);
  assert(e);
  center_ = Point{ival[0], ival[1], ival[2]};

  e = buffer.read(&scale_);
  assert(e);

  for (int i = 0; i < ct; ++i) {
    int idx{};
    size_t bts{};
    e = buffer.read((char *)&idx, sizeof(idx));
    assert(e);
    e = buffer.read((char *)&bts, sizeof(bts));
    assert(e);

    add_view(idx, bts, nullptr);
  }

  for (int i = 0; i < ct; ++i) {
    set_data(i, buffer.cursor());
    e = buffer.advance(view_bytes(i));
    assert(e);
  }
}


} // namespace dashmm

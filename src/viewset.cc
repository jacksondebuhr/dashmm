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


#include "dashmm/viewset.h"


/// \file
/// \brief Implementation of ViewSet


namespace dashmm {


void ViewSet::clear() {
  views_.clear();
  //n_digits_ = -1;
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

  retval += 3 * sizeof(int);  // for count and n_digits and role
  retval += 4 * sizeof(double); // for the center and the scale
  retval += views_.size() * sizeof(int);
  retval += views_.size() * sizeof(size_t);

  return retval;
}


void ViewSet::serialize(WriteBuffer &buffer) {
  bool e = buffer.write((char *)&role_, sizeof(role_));
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

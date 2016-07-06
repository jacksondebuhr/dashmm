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


#ifndef __DASHMM_VIEW_SET_H__
#define __DASHMM_VIEW_SET_H__


#include "dashmm/buffer.h"


namespace dashmm {


namespace {
  struct View {
    int index;
    size_t bytes;
    char *data;
  };
}


class ViewSet {
 public:
  ViewSet() : views_{} { }
  explicit ViewSet(int ndig) : views_{}, n_digits_{ndig} { }

  //destroy the represented data
  void destroy();

  //things to add to the set
  void add_view(int index);
  void add_view(int index, size_t bytes, char *data);

  //things to modify given entries
  void set_bytes(int view, size_t bytes) {views_[view].bytes = bytes;}
  void set_data(int view, char *data) {views_[view].data = data;}
  void set_n_digits(int ndig) {n_digits_ = ndig;}

  //things to read given entries
  int view_index(int view) const {return views_[view].index;}
  size_t view_bytes(int view) const {return views_[view].bytes;}
  char *view_data(int view) const {return views_[view].data;}
  int n_digits() const {return n_digits_;}

  //collective information
  int count() const {return views_.size();}
  size_t bytes() const; // this includes eveything read or written by the object

  //serialize and interpret - technically, deserialize is not
  // accurate. It will not copy the data, just collect where that data ought
  // to be. Hence the interpret name
  //
  // TODO: The error handling here is not very good. These fail assertions if
  // something goes wrong.
  void serialize(WriteBuffer &buffer);
  void interpret(ReadBuffer &buffer);

 private:
  std::vector<View> views_;
  int n_digits_;
};


} // namespace dashmm


#endif // __DASHMM_VIEW_SET_H__

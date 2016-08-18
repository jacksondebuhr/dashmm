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


/// \file include/dashmm/viewset.h
/// \brief Implementation of ViewSet


#include <vector>

#include "dashmm/buffer.h"
#include "dashmm/types.h"


namespace dashmm {


/// This type provides the basic data needed for a View
struct View {
  int index;
  size_t bytes;
  char *data;
};


/// A collection of Views of an Expansion
///
/// This object is intended to be used to set up subsets of the available
/// views in a given expansion, or to select a subset of the available views
/// for use in the DASHMM library.
class ViewSet {
 public:
  /// Create an empty ViewSet
  ViewSet()
    : views_{}, role_{kNoRoleNeeded}, center_{0.0, 0.0, 0.0}, scale_{1.0} { }

  /// Create an empty ViewSet, while setting some vitals
  ///
  /// \param role - the role of the represented views
  ViewSet(int ndig, ExpansionRole role, const Point &center, double scale)
    : views_{}, role_{role}, center_{center}, scale_{scale} { }

  /// Clear the ViewSet, including the n_digits and role.
  void clear();

  /// Add a view to the set by index.
  ///
  /// This will indicate that the view at the given index is part of the set,
  /// but the size of the view, and the location of the view data is unknown.
  /// This is used when querying an Expansion for a set of views, which will
  /// fill in the bytes and data members of the View.
  ///
  /// \param index - the index of the view to add
  void add_view(int index);

  /// Add full details of a view.
  ///
  /// This is used to add a complete description of a view to the ViewSet.
  ///
  /// \param index - the index of the view
  /// \param bytes - the size of the view's data
  /// \param data - the address of the data.
  void add_view(int index, size_t bytes, char *data);

  /// Set the data size of a given view
  void set_bytes(int view, size_t bytes) {views_[view].bytes = bytes;}

  /// Set the data pointer of a given view
  void set_data(int view, char *data) {views_[view].data = data;}

  /// Set the role of the view set
  void set_role(ExpansionRole role) {role_ = role;}

  /// Set the center of the view set
  void set_center(const Point &center) {center_ = center;}

  /// Set the scale for the expansion
  void set_scale(double s) {scale_ = s;}

  /// Get the index of a given view
  int view_index(int view) const {return views_[view].index;}

  /// Get the data size of a given view
  size_t view_bytes(int view) const {return views_[view].bytes;}

  /// Get the data pointer of a given view
  char *view_data(int view) const {return views_[view].data;}

  /// Get the role for the ViewSet
  ExpansionRole role() const {return role_;}

  /// Get the center for the ViewSet
  Point center() const {return center_;}

  /// Get the scale
  double scale() const {return scale_;}

  /// Get the number of views in this ViewSet
  int count() const {return views_.size();}

  /// Get the size of the data represented by this ViewSet
  ///
  /// This includes not only the sizes of the views themselves, but also the
  /// meta data used during serialization of that data.
  size_t bytes() const;

  /// Serialize the represented data into a buffer
  ///
  /// This is used predominantly for packing a set of Expansion Views into a
  /// parcel to send around the system by HPX-5. As such, this will likely
  /// never be used by DASHMM users.
  void serialize(WriteBuffer &buffer);

  /// Interpret from a buffer a ViewSet
  ///
  /// This will interpret the data in buffer as a ViewSet. This is the
  /// operation that is paired with serialize. The use case for this does
  /// not require a copy be made, and so this data is interpreted rather than
  /// read. See Buffer for the distinction.
  void interpret(ReadBuffer &buffer);

 private:
  std::vector<View> views_;
  ExpansionRole role_;
  Point center_;
  double scale_;
};


} // namespace dashmm


#endif // __DASHMM_VIEW_SET_H__

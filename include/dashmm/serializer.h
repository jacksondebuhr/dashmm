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

#ifndef __DASHMM_SERIALIZER_H__
#define __DASHMM_SERIALIZER_H__


/// \file
/// \brief Abstract interface for serializers


namespace dashmm {


class Serializer {
 public:
  /// Return the serial size of the given object
  ///
  /// \param object - pointer to the object in question
  ///
  /// \returns - the size in bytes of the serialized object
  virtual size_t size(void *object) const = 0;

  /// Serialize an object
  ///
  /// This will serialize the given object into the buffer provided. This will
  /// then return the address after the serialized data. So the return value
  /// will be buffer + size(object).
  ///
  /// \param object - the object to serialize
  /// \param buffer - the address in a buffer in which to serialize the object
  ///
  /// \returns - the next available byte in the buffer after having serialized
  ///            the given object.
  virtual void *serialize(void *object, void *buffer) const = 0;

  /// Deserialize an object
  ///
  /// This will deserialize the given object into the buffer provided. This
  /// will then return the address after the data used in the deserialization.
  /// The return value will, thus, be buffer + size(object).
  ///
  /// \param buffer - the buffer from which to deserialize the object
  /// \param object - the object that will receive the data
  ///
  /// \returns - the next unused byte in the buffer after having deserialized
  ///            the object at buffer.
  virtual void *deserialize(void *buffer, void *object) const = 0;
};


} // dashmm


#endif // __DASHMM_SERIALIZER_H__
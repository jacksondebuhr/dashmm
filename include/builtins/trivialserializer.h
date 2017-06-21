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

#ifndef __DASHMM_TRIVIAL_SERIALIZER_H__
#define __DASHMM_TRIVIAL_SERIALIZER_H__


/// \file
/// \brief Trivial Copy Serializer implementation


#include <cstring>

#include "dashmm/serializer.h"


namespace dashmm {


/// This is a serializer for trivial types
template <typename T>
class TrivialSerializer {
 public:
  size_t size(void *object) const override {
    return sizeof(T);
  }

  void *serialize(void *object, void *buffer) const override {
    memcpy(buffer, object, sizeof(T));
    return (reinterpret_cast<char *>(buffer) + sizeof(T));
  }

  void *deserialize(void *buffer, void *object) const override {
    memcpy(object, buffer, sizeof(T));
    return (reinterpret_cast<char *>(buffer) + sizeof(T));
  }
};


} // dashmm


#endif // __DASHMM_TRIVIAL_SERIALIZER_H__
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


#ifndef __DASHMM_BUFFER_H__
#define __DASHMM_BUFFER_H__


#include <cassert>
#include <cstring>


namespace dashmm {


class Buffer {
public:
  Buffer(char *base, size_t total)
    : base_{base}, offset_{base}, total_{total}, remain_{total} { }

  char *base() const {return base_;}
  size_t capacity() const {return total_;}

  char *cursor() const {return offset_;}
  size_t remain() const {return remain_;}

  bool complete() const {return remain_ == 0;}

  bool advance(size_t bytes) {
    if (complete()) return false;
    size_t adv = (bytes <= remain_) ? bytes : remain_;
    offset_ += adv;
    remain_ -= adv;
    return adv == bytes;
  }

protected:
  char *base_;
  char *offset_;
  size_t total_;
  size_t remain_;
};


class ReadBuffer : public Buffer {
 public:
  //
  ReadBuffer(char *base, size_t total) : Buffer{base, total} { }

  bool read(char *data, size_t bytes) {
    if (complete()) return false;

    size_t toread = (bytes <= remain_) ? bytes : remain_;
    memcpy(data, offset_, toread);

    offset_ += toread;
    remain_ -= toread;

    return toread == bytes;
  }

  template <typename T>
  bool read(T *value) {
    return read(reinterpret_cast<char *>(value), sizeof(T));
  }

  template <typename T>
  T *interpret() {
    T *retval = reinterpret_cast<T *>(cursor());
    if (!advance(sizeof(T))) {
      retval = nullptr;
    }
    return retval;
  }
};


class WriteBuffer : public Buffer {
 public:
  WriteBuffer(char *base, size_t total) : Buffer{base, total} { }

  bool write(const char *data, size_t bytes) {
    if (complete()) return false;

    size_t towrite = (bytes <= remain_) ? bytes : remain_;
    memcpy(offset_, data, towrite);

    offset_ += towrite;
    remain_ -= towrite;

    return towrite == bytes;
  }

  bool write(ReadBuffer &input) {
    //write either all of input, or until this is used up
    // return true only id all of input was written into this
    if (complete()) return false;

    size_t inputleft = input.remain();
    size_t towrite = (inputleft <= remain_) ? inputleft : remain_;
    bool e = input.read(offset_, towrite);
    assert(e);

    offset_ += towrite;
    remain_ -= towrite;

    return towrite == inputleft;
  }

  template <typename T>
  bool write(const T value) {
    return write(reinterpret_cast<const char *>(&value), sizeof(T));
  }
};


} // namespace dashmm


#endif // __DASHMM_BUFFER_H__

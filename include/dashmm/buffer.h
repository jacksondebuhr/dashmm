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


/// \file include/dashmm/buffer.h
/// \brief Utility classes to make serialization of Expansions easier


#include <cassert>
#include <cstring>


namespace dashmm {


/// The Buffer object represents a segment of memory.
///
/// The object takes in the base and size of the buffer, and will keep track
/// of a cursor in the buffer. There is little reason to use this class
/// by itself. Instead, the real utility is provided in ReadBuffer and
/// WriteBuffer.
class Buffer {
public:
  /// Create the buffer from a given address and size.
  Buffer(char *base, size_t total)
    : base_{base}, offset_{base}, total_{total}, remain_{total} { }

  /// Return the base of the buffer's memory
  char *base() const {return base_;}

  /// The total capacity of the buffer
  size_t capacity() const {return total_;}

  /// The current location of the cursor
  char *cursor() const {return offset_;}

  /// The number of bytes remaining in the buffer
  size_t remain() const {return remain_;}

  /// Has the buffer been exhausted?
  bool complete() const {return remain_ == 0;}

  /// Advances the cursor by bytes bytes.
  ///
  /// \param bytes - the number of bytes to advance the cursor
  ///
  /// \returns - true if the advance occurred completely; false otherwise,
  ///            including the case where the buffer was advanced some, but
  ///            not the full distance specified by the argument.
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


/// Buffer to read serialized data from
class ReadBuffer : public Buffer {
 public:
  /// Create a buffer for reading
  ReadBuffer(char *base, size_t total) : Buffer{base, total} { }

  /// Read data from the buffer
  ///
  /// This will read bytes bytes from the given data into the memory at
  /// the address data. Note that this makes a copy of the data. For an
  /// alternative, see interpret() below.
  ///
  /// \param data - the address into which the data will be copied
  /// \param bytes - the number of bytes to copy.
  ///
  /// \returns - true if bytes bytes are successfully read; false on error
  ///            or partial read.
  bool read(char *data, size_t bytes) {
    if (complete()) return false;

    size_t toread = (bytes <= remain_) ? bytes : remain_;
    memcpy(data, offset_, toread);

    offset_ += toread;
    remain_ -= toread;

    return toread == bytes;
  }

  /// Template overload of read for ease-of-use.
  template <typename T>
  bool read(T *value) {
    return read(reinterpret_cast<char *>(value), sizeof(T));
  }

  /// Interpret some data from the buffer as a given type
  ///
  /// This will return a pointer with template type into the buffer.
  /// This will advance the cursor of the buffer.
  ///
  /// NOTE: the compiler will not be able to perform type deduction for this
  /// routine, so the explicit type will need to be specified for uses of
  /// this method:
  ///
  /// int *code = buf.interpret<int>();
  ///
  /// \returns - pointer to the interpreted data; nullptr on error.
  template <typename T>
  T *interpret() {
    T *retval = reinterpret_cast<T *>(cursor());
    if (!advance(sizeof(T))) {
      retval = nullptr;
    }
    return retval;
  }
};


/// Buffer to serialize data into
class WriteBuffer : public Buffer {
 public:
  // Create a write buffer
  WriteBuffer(char *base, size_t total) : Buffer{base, total} { }

  /// Write bytes bytes of the given data into the buffer.
  ///
  /// This will make a copy of the data, and will advance the buffer's cursor.
  /// If there is insufficient remaining capacity for the given data, the
  /// data that does fit will be copied. However, in this case, this method
  /// will return false, indicating the error.
  ///
  /// \param data - the buffer from which to write data
  /// \param bytes - the number of bytes to write
  ///
  /// \returns - true on success; false otherwise
  bool write(const char *data, size_t bytes) {
    if (complete()) return false;

    size_t towrite = (bytes <= remain_) ? bytes : remain_;
    memcpy(offset_, data, towrite);

    offset_ += towrite;
    remain_ -= towrite;

    return towrite == bytes;
  }

  /// Write the contents of a ReadBuffer into this object
  ///
  /// This will read as many bytes from input as can be handles by the
  /// remaining capacity of this object. Only if the number of bytes written
  /// makes up the full remaining capacity of the ReadBuffer will this method
  /// return as successful.
  ///
  /// \param input - ReadBuffer to exhaust into this
  ///
  /// \returns - true on successful exhaustion of input; false otherwise
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

  /// Templated overload of write for ease-of-use
  template <typename T>
  bool write(const T value) {
    return write(reinterpret_cast<const char *>(&value), sizeof(T));
  }
};


} // namespace dashmm


#endif // __DASHMM_BUFFER_H__

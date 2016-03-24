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


#ifndef __DASHMM_SOURCE_REF_H__
#define __DASHMM_SOURCE_REF_H__


/// \file include/dashmm/sourceref.h
/// \brief SourceRef object definition


#include <cstdlib>

#include <hpx/hpx.h>


namespace dashmm {


/// Source concept in DASHMM
///
/// To qualify as a Source for DASHMM, a type should be trivially copyable,
/// and should provide at one members accessible by name:
///
/// Point position;
///
/// For specific expansions, there will be further requirements on the contents
/// of the Source type for it to be compatible with the given Expansion.
/// In particular, sources will always have some required notion of 'charge'
/// for the Expansion in question. Often, this will be as simple as
///
/// double charge;
///
/// but could be more complicated in other cases. For details, please see the
/// specific Expansions in question.
///
/// Failure to provide a type with the required members will result in
/// compilation errors.



/// Reference to a set of sources
///
/// This is a reference object, meaning that it refers to the Source data in
/// the GAS, but does not contain those data. As such, one can pass this
/// class by value without worry. Also, because this is a reference, this
/// object cannot be used to destroy the memory in the global address space to
/// which it refers.
///
/// This class is a template parameterized by the Source type.
template <typename Source>
class SourceRef {
 public:
  using source_t = Source;
  using sourceref_t = SourceRef<Source>;

  /// Default constructor.
  SourceRef()
      : data_{HPX_NULL}, n_{0}, n_tot_{0} { }

  /// Construct from a specific address and count.
  SourceRef(hpx_addr_t data, size_t n, size_t n_tot)
      : data_{data}, n_{n}, n_tot_{n_tot} { }

  /// Returns if the reference is valid
  bool valid() const {return data_ != HPX_NULL;}

  /// Get a reference to a slice of the current reference
  ///
  /// This will return a SourceRef to a consecutive chunk of the records that
  /// begin at an offset from the start of this reference and that will contain
  /// n entries. If the input arguments to this method are invalid, then an
  /// invalid reference will be returned.
  ///
  /// \param offset - the offset from the start of this reference
  /// \param n - the number of entries in the slice
  ///
  /// \returns - the resulting SourceRef; may be invalid.
  sourceref_t slice(size_t offset, size_t n) const {
    if (offset > n_) {
      return sourceref_t{HPX_NULL, 0, 0};
    }
    if (offset + n > n_) {
      return sourceref_t{HPX_NULL, 0, 0};
    }
    hpx_addr_t addr = hpx_addr_add(data_, sizeof(Source) * offset,
                                   sizeof(Source) * n_tot_);
    return sourceref_t{addr, n, n_tot_};
  }

  /// Returns the number of Source records referred to.
  size_t n() const {return n_;}

  /// Returns the total number of Source records in the underlying allocation.
  size_t n_tot() const {return n_tot_;}

  /// Returns the global address of the referred to data.
  hpx_addr_t data() const {return data_;}

 private:
  /// Address of the first Source referred to by this object
  hpx_addr_t data_;
  /// Number of Sources referred to by this object
  size_t n_;
  /// Total number of Sources in the underlying GAS allocation
  size_t n_tot_;
};


} // namespace dashmm


#endif // __DASHMM_SOURCE_REF_H__

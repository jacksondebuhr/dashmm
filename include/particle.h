#ifndef __DASHMM_PARTICLE_H__
#define __DASHMM_PARTICLE_H__


/// \file include/particle.h
/// \brief Source and Target particle types.


#include <complex>

#include <hpx/hpx.h>

#include "include/point.h"


namespace dashmm {

using dcomplex_t = std::complex<double>; 

/// The data needed for source particles.
struct Source {
  double charge;
  Point position;
};


/// Reference to a set of sources
///
/// This is a reference object, meaning that it refers to the Source data in
/// the GAS, but does not contain those data. As such, one can pass this
/// class by value without worry. The data referred to will be a single block
/// of GAS which contains a number of source records.
class SourceRef {
 public:
  /// Default constructor.
  SourceRef() : data_{HPX_NULL}, n_{0} { }

  /// Construct a reference from Source data.
  ///
  /// This will allocate GAS memory and copy the input sources into the
  /// GAS allocation, before setting this object up as a reference to that
  /// GAS memory.
  SourceRef(Source *sources, int n);

  /// Construct from a specific address and count.
  SourceRef(hpx_addr_t data, int n) : data_{data}, n_{n} { }


  /// Destroy the particle data in GAS.
  ///
  /// This is needed because this object is a reference, and the destruction of
  /// this oject only destroys the reference.
  void destroy();

  /// Returns the number of Source records referred to.
  int n() const {return n_;}

  /// Returns the global address of the referred to data.
  hpx_addr_t data() const {return data_;}

 private:
  hpx_addr_t data_;
  int n_;
};


/// The data needed for target locations.
struct Target {
  Point position;
  dcomplex_t phi;
  /// DASHMM rearranges the input locations so we store the index from the
  /// original data so that we can copy results into the correct place.
  size_t index;
};


/// Reference to target data
///
/// This object is a reference, but unlike SourceRef, it refers to an LCO.
/// The LCO is a user-defined LCO that handles the concurrent contribution to
/// the potential of the target points. The data of the LCO is the actual
/// Target data and a few additional pieces of management data.
class TargetRef {
 public:
  /// Construct a default object
  TargetRef() : data_{HPX_NULL}, n_{0} { }

  /// Construct a reference from an existing global address
  TargetRef(hpx_addr_t data, int n) : data_{data}, n_{n} { }

  /// Construct a reference from input targets. This will create the LCO
  TargetRef(Target *targets, int n);

  /// As a reference, destruction of the referred to object must be done
  /// via destroy.
  void destroy();

  /// The number of target points represented by the object
  int n() const {return n_;}

  /// The global address of the referred to object
  hpx_addr_t data() const {return data_;}

  /// Indicate to the underlying LCO that all operations have been scheduled
  void finalize() const;

  /// Indicate to the underlying LCO that it should expect \param num more
  /// operations
  void schedule(int num) const;

  /// Contribute a S->T operation to the referred targets
  ///
  /// \param type - the type of the expansion serving the S->T operation
  /// \param n - the number of sources
  /// \param sources - the sources themselves
  void contribute_S_to_T(int type, int n, Source *sources) const;

  /// Contribute a M->T operation to the referred targets
  ///
  /// \param type - the type of expansion serving the M->T operation
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  /// \param scale - scaling factor
  void contribute_M_to_T(int type, size_t bytes, void *data, 
                         double scale) const;

  /// Contribute a L->T operation to the referred targets
  ///
  /// \param type - the type of expansion serving the L->T operation
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  /// \param scale - scaling factor
  void contribute_L_to_T(int type, size_t bytes, void *data, 
                         double scale) const;

 private:
  hpx_addr_t data_;
  int n_;
};


} //namespace dashmm


#endif // __DASHMM_PARTICLE_H__

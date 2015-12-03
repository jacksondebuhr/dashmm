#ifndef __DASHMM_EXPANSION_REF_H__
#define __DASHMM_EXPANSION_REF_H__


/// \file include/expansionref.h
/// \brief Interface to Expansion reference object


#include <memory>
#include <vector>

#include <hpx/hpx.h>

#include "include/expansion.h"
#include "include/index.h"
#include "include/particle.h"
#include "include/point.h"


namespace dashmm {


/// Reference to an Expansion
///
/// This object is a reference to data in GAS. As such, it can be passed by
/// value without worry. The point of the object is to provide an interface
/// to the GAS data without having to worry about the HPX-5 side of things,
/// unless one is interested.
///
/// As a reference, it attempts to have as many of the same methods as the
/// underlying Expansion object. However, some of the function signatures are
/// different, and some are missing.
///
/// The referred to object is actually a user-defined LCO that manages the
/// potentially concurrent contribution to the expansion.
class ExpansionRef {
 public:
  /// Construct the expansion from a given global address.
  ExpansionRef(int type, hpx_addr_t addr) : type_{type}, data_{addr} { }

  /// Destroy the GAS data referred by the object.
  void destroy();

  /// Return the global address of the referred data.
  hpx_addr_t data() const {return data_;}

  /// Is the object currently referring to global data?
  bool valid() const {return data_ != HPX_NULL;}

  /// What type of expansion is this referring to.
  int type() const {return type_;}

  //TODO: these...
  // So all of these basically have the following pattern: once the LCO is
  // ready, an action will spawn to continue the data to a relevant other
  // action.
  //

  //NOTE: These do not need to wait on the expansion
  // we will have access to the SourceRef, which knows counts. SO here we
  // just get the source data, compute the S_to_M and then set the LCO
  //
  // TODO: This needs more thinking too...
  //
  // See next. This is generally called on the prototype, and it then
  // contributes to the local.
  std::unique_ptr<Expansion> S_to_M(Point center, SourceRef sources) const;
  std::unique_ptr<Expansion> S_to_L(Point center, SourceRef sources) const;

  //NOTE: These *do* have to wait for the expansion
  // This is a call when on the expansion containted, that will perform the
  // translation, and continue that with the correct code to the expansion LCO
  //
  //TODO: These need more thinking...
  //
  //The results of these are typically added to a particular target or
  // something.
  std::unique_ptr<Expansion> M_to_M(int from_child, double s_size) const;
  std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                    Index t_index) const;
  std::unique_ptr<Expansion> L_to_L(int to_child, double t_size) const;

  /// Apply the effect of a multipole expansion to targets
  ///
  /// This will compute the effect of this multipole expansion on the given
  /// @p targets. Note that this is an asynchronous operation. This will only
  /// perform work once this expansion is ready.
  ///
  /// \param targets - the target for which the multipole expansion is applied
  void M_to_T(TargetRef targets) const;

  /// Apply the effect of a local expansion to targets
  ///
  /// This will compute the effect of this local expansion on the given
  /// @p targets. Note that this is an asynchronous operation. This will
  /// only perform work once this expansion is ready.
  ///
  /// \param targets - the target for which the local expansion is applied
  void L_to_T(TargetRef targets) const;

  /// Apply effect of sources to targets
  ///
  /// This will compute the effect of the given @p sources on the given
  /// @p targets. Note that this is an asynchronous operation. This may
  /// return before the contribution to the targets has been computed.
  ///
  /// \param sources - a reference to the source points
  /// \param targets - a reference to the target points
  void S_to_T(SourceRef sources, TargetRef targets) const;

  /// Add the given expansion to this expansion
  ///
  /// This will add the @p summand to this expansion. This is an asynchronous
  /// operation, and will only complete once @p summand is set. This does
  /// schedule a contribution, so this should only be called before the call
  /// to finalize().
  ///
  /// \param summand - a reference to the expansion to add to this one
  void add_expansion(ExpansionRef summand);

  /// Create a new expansion of the same type referred to by this object.
  ///
  /// \param center - the center point for the next expansion.
  ///
  /// \returns - the resulting expansion.
  std::unique_ptr<Expansion> get_new_expansion(Point center) const;

  //TODO: methods to make the inputs easy
  // these will wrap up the HPX stuff so the user can do "obvious" seeming
  // thigns instead.

  //TODO: Should this be private? We don't expect users to use this...
  /// Contribute to the referred expansion
  ///
  /// This will setup the given @p payload with the correct internal code
  /// and will call the appropriate set operation on the referred LCO. This
  /// will result in the add_expansion method of the expansion being called.
  ///
  /// \param bytes - the size of the input serialized expansion
  /// \param payload - the serialized expansion data
  void contribute(size_t bytes, char *payload);

  /// Signal to the expansion that all operations have been scheduled.
  ///
  /// This should be called only after all possible contributions to the
  /// expansion have been scheduled. The underlying LCO cannot trigger until
  /// finalize() has been called, and so any work dependent on this expansion
  /// would be on-hold until the call to finalize().
  void finalize() const;

  /// Signal to the expansion that it should expect an operation
  ///
  /// This will inform the underlying LCO of another eventual contribution to
  /// its value. The asynchronous nature of the computation in DASHMM means that
  /// contributions will happen when they are ready. schedule() will inform
  /// that a contribution is on the way.
  void schedule() const;

 private:
  int type_;
  hpx_addr_t data_;     //this is the LCO
};


/// Create a global version of a given expansion
///
/// This will export the local data represented by the given expansion into
/// the global address space, and return an ExpansionRef to that data. This
/// mechanism relies on the release() method of the particular expansion to
/// provide a packed representation of the relevant data for the expansion.
/// For more details, see the DASHMM Advanced User Guide.
ExpansionRef globalize_expansion(std::unique_ptr<Expansion> exp);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__

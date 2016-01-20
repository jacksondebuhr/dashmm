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
#include "include/types.h"


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
  ExpansionRef(int type, hpx_addr_t addr, int n_digits) 
    : type_{type}, data_{addr}, n_digits_{n_digits} { }

  /// Destroy the GAS data referred by the object.
  void destroy();

  /// Return the global address of the referred data.
  hpx_addr_t data() const {return data_;}

  /// Is the object currently referring to global data?
  bool valid() const {return data_ != HPX_NULL;}

  /// What type of expansion is this referring to.
  int type() const {return type_;}

  /// Accuracy of expansion
  int accuracy() const {return n_digits_;}

  /// Set this expansion with the multipole expansion of the given sources
  ///
  /// This will set the expansion with the multipole expansion computed for
  /// the given @p sources, and with the given @p center. Note that this is
  /// an asynchronous operation. The contribution will be scheduled, but will
  /// not necessarily be complete when this function returns.
  ///
  /// Note that this does not set the expansion to be equal to the computed
  /// multipole expansion. Instead, it adds the computed multipole to the
  /// current contents of the expansion. Further, this must not be called
  /// after finalize().
  ///
  /// \param center - the center of the computed expansion
  /// \param sources - a reference to the sources from which to compute the
  ///                  multipole expansion
  /// \param scale - scaling factor
  void S_to_M(Point center, SourceRef sources, double scale) const;

  /// Set this expansion with the local expansion of the given sources
  ///
  /// This will set the expansion with the local expansion computed for
  /// the given @p sources, and with the given @p center. Note that this is
  /// an asynchronous operation. The contribution will be scheduled, but will
  /// not necessarily be complete when this function returns.
  ///
  /// Note that this does not set the expansion to be equal to the computed
  /// local expansion. Instead, it adds the computed local to the
  /// current contents of the expansion. Further, this must not be called
  /// after finalize().
  ///
  /// \param center - the center of the computed expansion
  /// \param sources - a reference to the sources from which to compute the
  ///                  local expansion
  /// \param scale - scaling factor
  void S_to_L(Point center, SourceRef sources, double scale) const;

  /// Contribute a translated multipole moment to this expansion
  ///
  /// This will translate the given multipole expansion and then add it to
  /// this expansion. This is an asynchronous operation; this can return
  /// before the contribution to this expansion has been made. This must not
  /// be called after finalize().
  ///
  /// \param source - the multipole expansion to translate
  /// \param from_child - the child from which the expansion will occur
  /// \param s_size - the size of the child node for @p source
  void M_to_M(ExpansionRef source, int from_child, double s_size) const;

  /// Contribute a translated multipole moment to this local expansion
  ///
  /// This will translate a multipole expansion into a local expansion and
  /// add it to this expansion. This operation is asynchronous; this can
  /// return before the contribution has been made. This must not be called
  /// after finalize().
  ///
  /// \param source - the multipole expansion to translate
  /// \param s_index - the index of the node containing @p source
  /// \param s_size - the size of the node containing @p source
  /// \param t_index - the index of the node containing this expansion
  void M_to_L(ExpansionRef source, Index s_index, double s_size,
              Index t_index) const;

  /// Contribute a translated local expansion to this local expansion
  ///
  /// This will translate the given local expansion into a local expansion
  /// that can be added to this expansion. This operation is asynchronous;
  /// this can return before the contribution to this expansion has been made.
  /// This must not be called after finalize().
  ///
  /// \param source - the local expansion to translate
  /// \param to_child - the child to which the expansion is being translated
  /// \param t_size - the size of the child node
  void L_to_L(ExpansionRef source, int to_child, double t_size) const;

  /// Apply the effect of a multipole expansion to targets
  ///
  /// This will compute the effect of this multipole expansion on the given
  /// @p targets. Note that this is an asynchronous operation. This will only
  /// perform work once this expansion is ready.
  ///
  /// \param targets - the target for which the multipole expansion is applied
  /// \param scale - scaling factor
  void M_to_T(TargetRef targets, double scale) const;

  /// Apply the effect of a local expansion to targets
  ///
  /// This will compute the effect of this local expansion on the given
  /// @p targets. Note that this is an asynchronous operation. This will
  /// only perform work once this expansion is ready.
  ///
  /// \param targets - the target for which the local expansion is applied
  /// \param scale - scaling factor
  void L_to_T(TargetRef targets, double scale) const;

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
  /// to finalize(). This routine takes ownership of the supplied expansion,
  /// and will free any resources associated with the expansion.
  ///
  /// \param summand - a reference to the expansion to add to this one
  void add_expansion(ExpansionRef summand);

  /// Create a new expansion of the same type referred to by this object.
  ///
  /// \param center - the center point for the next expansion.
  ///
  /// \returns - the resulting expansion.
  std::unique_ptr<Expansion> get_new_expansion(Point center, 
                                               int n_digits) const;

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

  /// Contribute to the referred expansion
  ///
  /// This will setup the given @p payload with the correct internal code
  /// and will call the appropriate set operation on the referred LCO. This
  /// will result in the add_expansion method of the expansion being called.
  ///
  /// \param bytes - the size of the input serialized expansion
  /// \param payload - the serialized expansion data
  void contribute(size_t bytes, char *payload);

 private:
  int type_;
  hpx_addr_t data_;     //this is the LCO
  int n_digits_; 
};


/// Create a global version of a given expansion
///
/// This will export the local data represented by the given expansion into
/// the global address space, and return an ExpansionRef to that data. This
/// mechanism relies on the release() method of the particular expansion to
/// provide a packed representation of the relevant data for the expansion.
/// For more details, see the DASHMM Advanced User Guide.
///
/// \param exp - the expansion to globalize
/// \param where - an HPX address indicating an address which should be local
///                to the globalized expansion
///
/// \returns - a reference to the resulting expansion
ExpansionRef globalize_expansion(std::unique_ptr<Expansion> exp,
                                 hpx_addr_t where);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__

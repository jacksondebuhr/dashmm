#ifndef __DASHMM_EXPANSION_H__
#define __DASHMM_EXPANSION_H__


/// \file include/expansion.h
/// \brief Abstract interface to expansions


#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include <hpx/hpx.h>

#include "include/index.h"
#include "include/particle.h"
#include "include/types.h"


namespace dashmm {

using dcomplex_t = std::complex<double>; 

/// The abstract interface for expansions usable in DASHMM
///
/// This interface specifies the requirements for expansions that a user
/// may add to DASHMM. In general, the user will need to implement the
/// expansion and then somewhere in their program call register_expansion()
/// providing the creation and interpretation functions that the user also
/// implements. For details, see the DASHMM Advanced User Guide.
class Expansion {
 public:
  virtual ~Expansion() { }

  /// Release the internal data for an expansion
  ///
  /// Expansion objects need to support the ability to export the data making
  /// up the expansion in a contiguous chunk of memory. These serialized
  /// versions of the expansion are used throughout the system. When release
  /// is called, the expansion object is then 'empty' in a sense, and the
  /// expansion no longer owns any data.
  ///
  /// At the beginning of the serialized version of the expansion must be
  /// two integers. The first is a code that will be used internally by
  /// DASHMM and does not need to be given a value. The second is the integer
  /// type of this expansion, which should be given a value by the user.
  ///
  /// The caller of this function assumes ownership of the released data,
  /// and it is the caller's responsibility to free() that data. Note,
  /// this means the data should be allocated with malloc().
  virtual void *release() = 0;

  /// The size in bytes of the serialized expansion data
  virtual size_t bytes() const = 0;

  /// Returns if the expansion is valid.
  ///
  /// An expansion is valid if it has data associated with it. After calling
  /// release(), an expansion will be invalidated.
  virtual bool valid() const = 0;

  /// Returns the type identifier of the expansion
  ///
  /// Each expansion must have a unique type in the range
  /// [kFirstUserExpansionType, kLastUserExpansionType]. The user is responsible
  /// for managing the uniquess of user-defined types. For details see the
  /// DASHMM Advanced User Guide.
  virtual int type() const = 0;

  /// Returns the accuracy of the expansion
  virtual int accuracy() const = 0; 

  /// Indicates if the expansion provides local expansion operators.
  virtual bool provides_L() const = 0;

  /// Indicates if the expansion provides exponential expansion operators.
  virtual bool provides_exp() const = 0;

  /// The number of terms in the expansion.
  virtual size_t size() const = 0;
  
  /// The point around which the expansion is defined.
  virtual Point center() const = 0;

  /// Get a term of the expansion.
  ///
  /// The input should be in the range [0, size()). The particular
  /// implementation will decide if this is range checked. To be safe, assume
  /// that the input is not range checked.
  ///
  /// \param i - the term to return
  ///
  /// \returns - the complex double version of the term; real expansions should
  ///            nevertheless return terms as a complex number for compatibility
  //             with complex expansions.
  virtual dcomplex_t term(size_t i) const = 0;

  /// Create a multipole expansion for a given set of source points
  ///
  /// This uses the given sources to create a multipole expansion centered
  /// at \p center. The sources are provided as pointers to the \p first and
  /// one past the \p last source record.
  ///
  /// \param center - the point around which to form the expansion
  /// \param first - address of the first source
  /// \param last - address of one past the last source
  /// \param scale - scaling factor
  ///
  /// \returns - The resulting multipole expansion
  virtual std::unique_ptr<Expansion> S_to_M(Point center, 
                                            Source *first, Source *last, 
                                            double scale) const = 0;

  /// Create a local expansion for a given set of source points
  ///
  /// This uses the given sources to create a local expansion centered
  /// at \p center. The sources are provided as pointers to the \p first and
  /// one past the \p last source record.
  ///
  /// \param center - the point around which to form the expansion
  /// \param first - address of the first source
  /// \param last - address of one past the last source
  /// \param scale - scaling factor
  ///
  /// \returns - The resulting local expansion
  virtual std::unique_ptr<Expansion> S_to_L(Point center, 
                                            Source *first, Source *last,
                                            double scale) const = 0;

  /// Change center of a multipole expansion
  ///
  /// This will convert a multipole expansion of a child node to a multipole
  /// expansion for a parent node. The change in center is specified by
  /// indicating from which child the expansion is being shifted.
  ///
  /// \param from_child - the child index from which the expansion is being
  ///                     converted.
  /// \param s_size - the size of the child node
  ///
  /// \returns - The resulting multipole expansion
  virtual std::unique_ptr<Expansion> M_to_M(int from_child,
                                            double s_size) const = 0;

  /// Convert a multipole to a local expansion
  ///
  /// This will convert a multipole expansion in the source tree into
  /// a local expansion in the target tree.
  ///
  /// \param s_index - the index specifying the source node
  /// \param s_size - the size of the source node
  /// \param t_index - the index specifying the target node
  ///
  /// \returns - The resulting local expansion.
  virtual std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                             Index t_index) const = 0;

  /// Convert a local expansion to a local expansion for a child
  ///
  /// This converts the local expansion of a parent node into a local expansion
  /// valid for a child of the parent.
  ///
  /// \param to_child - the index of the child to which the expansion is being
  ///                   transformed.
  /// \param t_size - the size of the child node.
  ///
  /// \returns - The resulting local expansion.
  virtual std::unique_ptr<Expansion> L_to_L(int to_child,
                                            double t_size) const = 0;

  /// Apply a multipole expansion to a set of targets
  ///
  /// \param first - the first target point
  /// \param last - one past the last target point
  /// \param scale - scaling factor
  virtual void M_to_T(Target *first, Target *last, double scale) const = 0;

  /// Apply a local expansion to a set of targets
  ///
  /// \param first - the first target point
  /// \param last - one past the last target point
  /// \param scale - scaling factor
  virtual void L_to_T(Target *first, Target *last, double scale) const = 0;

  /// Compute the direct interaction between sources and targets
  ///
  /// \param s_first - the first source point
  /// \param s_last - one past the last source point
  /// \param t_first - the first target point
  /// \param t_last - one past the last target point
  virtual void S_to_T(Source *s_first, Source *s_last,
                      Target *t_first, Target *t_last) const = 0;

  /// Add an expansion to this expansion
  ///
  /// Given another expansion (assumed to be of the same type), add the
  /// expansions together.
  virtual void add_expansion(const Expansion *temp1) = 0;

  /// Generate a new expansion of the same type
  ///
  /// This factory method allows an unspecified Expansion produce another of
  /// the correct type.
  ///
  /// \param center - the center point of the expansion
  ///
  /// \returns - the new expansion
  virtual std::unique_ptr<Expansion> get_new_expansion(Point center) const = 0;
};


/// The signature of a function that creates an expansion
///
/// To register a user-defined expansion with DASHMM, the user will need to
/// provide a function that can be used to create that sort of expansion
/// given a Point indicating the center of the expansion.
typedef Expansion *(*expansion_creation_function_t)(double, double, double, 
                                                    int);

/// The signature of a function that interprets existing data as an expansion
///
/// To register a user-defined expansion with DASHMM, the user will need to
/// provide a function that can be used to interpret exisiting data as an
/// expansion of the new type.
///
/// Some expansion methods do not need access to the particular data for the
/// expansion, so this function must be able to support creating a 'thin' or
/// 'empty' expansion of the new type. One example of this is the S_to_T
/// function. The behavior of this function is the same for every instance of
/// the same type of expansion.
typedef Expansion *(*expansion_interpret_function_t)(void *data, size_t bytes, 
                                                     int n_digits); 


/// Register a user-defined expansion with DASHMM
///
/// This will connect the two specified functions to the given \p type for
/// DASHMM. This will allow the system to interact with user-defined expansions.
/// For details, see the DASHMM Advanced User Guide.
///
/// \param type - the type identifier of the expansion being registered.
/// \param creator - the action identifier of the creation function for the
///                 type being registered.
/// \param interpreter - the action identifier of the interpretation function
///                      for the type being registered.
///
/// \returns - kSuccess on success; kDomainError otherwise indicating that the
///            specified type is either out of the user range, or already taken.
ReturnCode register_user_expansion(int type, hpx_action_t creator,
                                  hpx_action_t interpreter);


//NOTE: not intended for end-user use


/// Register a user-defined expansion with DASHMM
///
/// This will connect the two specified functions to the given \p type for
/// DASHMM. This will allow the system to interact with user-defined expansions.
/// For details, see the DASHMM Advanced User Guide.
///
/// \param type - the type identifier of the expansion being registered.
/// \param creator - the action identifier of the creation function for the
///                 type being registered.
/// \param interpreter - the action identifier of the interpretation function
///                      for the type being registered.
/// \param user - nonzero if the expansion is a user defined expansion; 0
///               for built-in expansions
///
/// \returns - kSuccess on success; kDomainError otherwise indicating that the
///            specified type is either out of the user range, or already taken.
ReturnCode register_expansion(int type, hpx_action_t creator,
                                  hpx_action_t interpreter, int user);


/// Initialize the expansion registration table
void init_expansion_table();


/// Finalize (clean up) the expansion registration table
void fini_expansion_table();


/// Interpret some data as an expansion object
///
/// This will interpret the given data as an expansion of the given type.
/// The returned expansion will be constructed from the interpret function
/// provided by the user during expansion registration.
///
/// Note that the resulting expansion will assume ownership of the data
/// passed in. When the returned object is destroyed, it will attempt to
/// free the data. If this is an error, the user will need to call release()
/// before letting the object be destroyed.
///
/// Note also that the data might be a nullptr. This is to allow for shallow
/// construction of an Expansion when the specific data for the expansion is
/// not needed. Examples of this would be to get the size(), or to perform
/// S_to_T().
///
/// \param type - the type identifier for the expansion
/// \param data - the data for the expansion, may be nullptr
/// \param size - the size of the data for the expansion
///
/// \returns - the resulting expansion
std::unique_ptr<Expansion> interpret_expansion(int type, void *data,
                                               size_t size, int n_digits);

/// Create an expansion of the given type
///
/// This will create a new expansion of the given type. The returned expansion
/// will be constructed from the create function provided by the user during
/// expansion registration.
///
/// \param type - the type identifier for the expansion
/// \param center - the Point around which to create the expansion
///
/// \returns - the resulting expansion
std::unique_ptr<Expansion> create_expansion(int type, Point center, 
                                            int n_digits); 


} // namespace dashmm


#endif // __DASHMM_EXPANSION_H__

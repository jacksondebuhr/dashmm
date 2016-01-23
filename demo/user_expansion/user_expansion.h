// The contents of this file (and of user_expansion.cc) are intended to be a
// skeleton of what is required to define an Expansion subclass that can be
// registered with DASHMM. The comments in these files will explain what is
// needed and where the user will have to implement their specific case.


#ifndef __USER_EXPANSION_H__
#define __USER_EXPANSION_H__


#include <memory>

#include "dashmm.h"


// Forward declaration of the User Expansion's data. It is declared in full
// in user_expansion.cc.
struct UserData;


// To register an expansion with DASHMM, it must have an integral identifier.
// DASHMM  requires that this type identifier falls in the range between
// dashmm::kFirstUserExpansionType and dashmm::kLastUserExpansionType
// inclusive. The current values of these constants allow for 1000 user
// registered Expansion subclasses.
//
// If the user is registering multiple Expansion subclasses they
// are responsible for assuring a unique identifier. Of course, DASHMM will
// complain if an attempt is made to register the same identifier twice.
constexpr int kUserExpansionType = dashmm::kFirstUserExpansionType;


// A user-defined expansion will be a subclass of dashmm::Expansion, an
// abstract base class defining the required interface.
//
// When specializing this to a specific need, there is likely a better name
// than 'User', and so that should be changed to whatever is appropriate.
//
// The implementation strategy of Expansion objects is similar to a handle to
// some other memory. The reason for this is that DASHMM will need to send
// serialized versions of the expansion data across the system, and so
// Expansion objects are required to have the capability of releasing their
// data in a contiguous block of memory on request (see release() below).
class User : public dashmm::Expansion {
 public:
  // Inside the library there are a few ways that Expansion objects are
  // constructed. The first is from the expansion point and the accuracy
  // requirement (this is called n_digits here, but it can be whatever number
  // is most convenient for the specific case). This constructor will create
  // a new allocation of memory for the UserData object that this object will
  // handle.
  User(dashmm::Point center, int n_digits);

  // The second creation method is via capturing some existing UserData and
  // causing this object to handle that memory. As the serialized Expansion
  // data is moved around the system, it will be necessary to interpret that
  // data as the correct subclass of Expansion. This constructor is where this
  // occurs.
  //
  // There are some operations that do not need the actual expansion terms
  // (such as S->T), and so this constructor should also work for
  // ptr == nullptr. The object so constructed will not own any data, but will
  // nevertheless be able to perform those operations not requiring the
  // expansion data. In these 'empty' constructions, n_digits will still be
  // supplied.
  User(UserData *ptr, size_t bytes, int n_digits)
      : data_{ptr}, bytes_{bytes}, acc_{n_digits} { }

  // Expansion objects that still own their data should release that resource
  // when they are destroyed. See release().
  ~User();

  // This routine will supply the serialized version of the data in this
  // expansion object. After a call to release(), this object no longer owns
  // the associated data, and will not free it when this object is destroyed.
  //
  // That this interface is required is the reason that these objects should
  // be implemented as a handle. Otherwise, this method would need to first
  // copy the data into the serialized form before releasing the data. The
  // serialized versions of Expansion objects are interpreted frequently, and
  // are released frequently. So performing these extra copies would impact
  // performance. That being said, one can implement an Expansion in other
  // style, but one should be aware that there could be performance impacts.
  void *release() override;

  // The size of the serialized expansion. This is especially important for
  // Expansions that might have different numbers of terms based on the
  // supplied accuracy parameter.
  size_t bytes() const override {return bytes_;}

  // This tests if this object is valid; that is, does this object own some
  // data.
  bool valid() const override {return data_ != nullptr;}

  // This is the type identifier for this Expansion subclass.
  int type() const override {return kUserExpansionType;}

  // This returns the accuracy paramter the object was constructed with.
  int accuracy() const override {return acc_;}

  // This indicates if this Expansion provides local translations of the
  // multipole expansions. To use FMM, the Expansion must provide L operations.
  bool provides_L() const override {return true;}

  // This gives the number of terms in the expansion.
  size_t size() const override;

  // This gives the point around which this Expansion is expanded.
  dashmm::Point center() const override;

  // This returns a specific term in the expansion. dcomplex_t is an
  // alias for std::complex<double>.
  dcomplex_t term(size_t i) const override;

  // This routine will generate the multipole expansion about the given
  // center for the given points.
  std::unique_ptr<dashmm::Expansion> S_to_M(dashmm::Point center,
                                            dashmm::Source *first,
                                            dashmm::Source *last,
                                            double scale) const override;

  // This will generate a local expansion at the given center for the
  // given points.
  std::unique_ptr<dashmm::Expansion> S_to_L(dashmm::Point center,
                                            dashmm::Source *first,
                                            dashmm::Source *last,
                                            double scale) const override;

  // This will translate a multipole expansion from a child node to its parent.
  std::unique_ptr<dashmm::Expansion> M_to_M(int from_child,
                                            double s_size) const override;

  // This will translate a multipole expansion into a local expansion.
  std::unique_ptr<dashmm::Expansion> M_to_L(dashmm::Index s_index,
                                        double s_size,
                                        dashmm::Index t_index) const override;

  // This will translate a local expansion from a parent node to one of
  // its children.
  std::unique_ptr<dashmm::Expansion> L_to_L(int to_child,
                                            double t_size) const override;

  // This will compute the effect of a multipole expansion on a set of
  // target points.
  void M_to_T(dashmm::Target *first, dashmm::Target *last,
              double scale) const override;

  // This will compute the effect of a local expansion on a set of target
  // points.
  void L_to_T(dashmm::Target *first, dashmm::Target *last,
              double scale) const override;

  // This will compute the effect of a set of sources on a set of target points.
  void S_to_T(dashmm::Source *s_first, dashmm::Source *s_last,
              dashmm::Target *t_first, dashmm::Target *t_last) const override;

  // The routine will add the provided expansion to this expansion.
  void add_expansion(const dashmm::Expansion *temp1) override;

  // There are a few cases where DASHMM needs to create a new object of the
  // subclass being used without knowing the type of that subclass (it has
  // a pointer to the base class). This method will create a new object of
  // the correct type.
  std::unique_ptr<dashmm::Expansion> get_new_expansion(
                                          dashmm::Point center) const override;

 private:
  // This object acts as a handle to some other piece of memory. The details
  // of the UserData type can be found in user_expansion.cc.
  UserData *data_;

  // This stores the size (in bytes) of data_.
  size_t bytes_;

  // This stores the accuracy parameter passed into the constuctor for this
  // object.
  int acc_;
};


// To register an Expansion with DASHMM, one has to supply a function that
// will create a new Expansion of the registered type. Generally, this
// function should construct the object via the first constructor for this
// Expansion type. Note that this will return the new object as a pointer
// to the base class.
dashmm::Expansion *create_user_expansion(double cx, double cy, double cz,
                                         int acc);

// To register an Expansion with DASHMM, one has to supply a function that can
// interpret an existing serialization of an Expansion as an object of the
// correct type. Generally, this function should construct the object via the
// second constructor for this Expansion type. Note that this will return the
// new object as a pointer to the base class.
dashmm::Expansion *interpret_user_expansion(void *data, size_t bytes, int acc);


#endif // __USER_EXPANSION_H__

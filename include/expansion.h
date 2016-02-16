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



/// To qualify for the Expansion concept in DASHMM, a class must satisfy the
/// following criteria. This file is not included anywhere in DASHMM, but it
/// is in the source distribution as an example, and to explain the
/// Expansion concept.


/// Expansions in DASHMM are template classes parameterized over the types of
/// the sources and targets. A full description of the requirements of the
/// Source and Target types can be found elsewhere.
///
/// When creating a user-defined Expansion, the name Expansion in the following
/// should be replaced by the name of the new Expansion type.

template <typename Source, typename Target>
class Expansion {
 public:
  /// Expansions are required to alias their serialized data type as
  /// 'contents_t' so that internally DASHMM can cast message buffers into the
  /// correct type as soon as possible.
  using contents_t = ExpansionData;

  //TODO: decide if this is actually useful
  using source_t = Source;
  using target_t = Target;

  /// Expansions must provide two constructors
  ///
  /// The first creates the object with the given center and the given
  /// accuracy parameter. For FMM, this might be the number of digits
  /// requested. The type is not obligated to use either of these parameters,
  /// but the type must provide a constructor of this form.
  Expansion(Point center, int n_digits);

  /// The second creates the expansion from previously existing data. The
  /// bytes argument allows for variable length expansions in DASHMM.
  ///
  /// Further, this constructor needs to be able to operate in a 'shallow'
  /// mode, where ptr = nullptr. This allows for situations where Expansion
  /// is needed, but the speicif data for the expansion is not. The exemplar of
  /// this use is to perform an S->T operation.
  ///
  /// Again, n_digits is required, but could be ignored.
  ///
  Expansion(contents_t *ptr, size_t bytes, int n_digits)

  /// The destructor should delete the allocated memory of the object. In the
  /// simplest style of implementation, this means that only if the object
  /// is valid() will the destructor delete any memory.
  ~Expansion();

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
  /// given as the n_digits in the construction of the specific instance of
  /// the object.
  ///
  /// The caller of this function assumes ownership of the released data,
  /// and it is the caller's responsibility to release the data. This should be
  /// done with delete [], as the memory will have been allocated with
  /// new char [some_size].
  ///
  /// After having been released, this object must report false when valid() is
  /// invoked.
  ///
  /// The simplest way to implement this is to have the object store a pointer
  /// to memory allocated on the heap (cast from new char [size]). And then
  /// release() will just return that pointer, and set this objects member to
  /// nullptr.
  void *release();

  /// The size in bytes of the serialized expansion data
  size_t bytes() const;

  /// Returns if the expansion is valid.
  ///
  /// An expansion is valid if it has data associated with it. After calling
  /// release(), an expansion will be invalidated.
  bool valid() const;

  /// Returns the accuracy of the expansion - this is the number given to the
  /// constructor as n_digits, or is some nonsense value if n_digits is not
  /// used.
  int accuracy() const;

  /// The number of terms in the expansion.
  size_t size() const;

  /// The point around which the expansion is defined.
  Point center() const;

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
  ///            with complex expansions.
  dcomplex_t term(size_t i) const;


  /// In the following operators, it is not required that all are implemented.
  /// Only those being used by the Method are needed. In practice, all of the
  /// operators with an 'M' in their title are always implemented. And only for
  /// FMM-like methods are the 'L' operations needed. Of course, S->T is always
  /// required.


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
  void S_to_M(Point center, Source *first, Source *last,
              double scale) const;

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
  std::unique_ptr<Expansion> S_to_L(Point center,
                                    Source *first, Source *last,
                                    double scale) const;

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
  std::unique_ptr<Expansion> M_to_M(int from_child,
                                    double s_size) const;

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
  std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                    Index t_index) const;

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
  std::unique_ptr<Expansion> L_to_L(int to_child,
                                    double t_size) const;

  /// Apply a multipole expansion to a set of targets
  ///
  /// \param first - the first target point
  /// \param last - one past the last target point
  /// \param scale - scaling factor
  void M_to_T(Target *first, Target *last, double scale) const;

  /// Apply a local expansion to a set of targets
  ///
  /// \param first - the first target point
  /// \param last - one past the last target point
  /// \param scale - scaling factor
  void L_to_T(Target *first, Target *last, double scale) const;

  /// Compute the direct interaction between sources and targets
  ///
  /// \param s_first - the first source point
  /// \param s_last - one past the last source point
  /// \param t_first - the first target point
  /// \param t_last - one past the last target point
  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const;

  /// Add an expansion to this expansion
  ///
  /// Given another expansion (assumed to be of the same type), add the
  /// expansions together.
  ///
  void add_expansion(const Expansion *temp1);
};

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
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;

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
  /// mode, where views contains no views. This allows for situations where
  /// Expansion is needed, but the specific data for the expansion is not.
  /// The exemplar of this use is to perform an S->T operation. The expansion
  /// will get a value for n_digits from the provided ViewSet.
  Expansion(ViewSet &views);

  /// The destructor should delete the allocated memory of the object. In the
  /// simplest style of implementation, this means that only if the object
  /// is valid() for any of its views should the object delete any data.
  ~Expansion();

  /// Release the internal data for an expansion
  ///
  /// Expansion objects need to support the ability to export the data making
  /// up the expansion into a given buffer (see view() below). Further,
  /// Expansion objects need to be able to be constructed in a shallow way
  /// from existing data. release() will break the association of the object
  /// with the data serving one or more valid views. When an Expansion is
  /// created from existing data, that data is owned by some other object or
  /// context, and for this object to free that memory would be an error.
  ///
  /// After having been released, this object must report false when valid() is
  /// invoked.
  ///
  /// The simplest way to implement this is to have the object store a pointer
  /// to memory allocated on the heap (cast from new char [size]). And then
  /// release() will just set this objects member to nullptr.
  void release();

  /// Returns if the indicated views are valid
  ///
  /// An expansion is valid if it has data associated with it. After calling
  /// release(), an expansion will be invalidated.
  bool valid(const ViewSet &view) const;

  /// Return the current number of views for the object.
  ///
  /// This will either be the full number for a directly created object, or
  /// it will be a smaller number if the object is created from a ViewSet.
  int view_count() const;

  /// Fill in the data for a view object
  ///
  /// Given a ViewSet that contains the view indices, this will populate the
  /// bytes and data members for the ViewSet. Note that this does not mean a
  /// copy will be made of the data, but merely that the data pointers will be
  /// valid.
  void get_views(ViewSet &view) const;

  /// Get all the current views of the object
  ViewSet get_all_views() const;

  /// Returns the accuracy of the expansion - this is the number given to the
  /// constructor as n_digits, or is some nonsense value if n_digits is not
  /// used.
  int accuracy() const;

  /// The point around which the expansion is defined.
  Point center() const;

  /// The number of terms in the expansion.
  size_t view_size(int view) const;

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
  dcomplex_t view_term(int view, size_t i) const;


  /// Create a multipole expansion for a given set of source points
  ///
  /// This uses the given sources to create a multipole expansion centered
  /// at \p center. The sources are provided as pointers to the \p first and
  /// one past the \p last source record.
  ///
  /// \param center - the point around which to form the expansion
  /// \param first - address of the first source
  /// \param last - address of one past the last source
  /// \param scale - scaling factor (one over the size of the source box)
  std::unique_ptr<expansion_t> S_to_M(Point center, source_t *first,
                                      source_t *last, double scale) const;

  /// Create a local expansion for a given set of source points
  ///
  /// This uses the given sources to create a local expansion centered
  /// at \p center. The sources are provided as pointers to the \p first and
  /// one past the \p last source record.
  ///
  /// \param center - the point around which to form the expansion
  /// \param first - address of the first source
  /// \param last - address of one past the last source
  /// \param scale - scaling factor (one over the size of the target box)
  ///
  /// \returns - The resulting local expansion
  std::unique_ptr<expansion_t> S_to_L(Point center,
                                      source_t *first, source_t *last,
                                      double scale) const;

  // TODO
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
  std::unique_ptr<expansion_t> M_to_M(int from_child,
                                      double s_size) const;

  // TODO
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
  std::unique_ptr<expansion_t> M_to_L(Index s_index, double s_size,
                                      Index t_index) const;

  // TODO
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
  std::unique_ptr<expansion_t> L_to_L(int to_child,
                                      double t_size) const;

  // TODO
  /// Apply a multipole expansion to a set of targets
  ///
  /// \param first - the first target point
  /// \param last - one past the last target point
  /// \param scale - scaling factor (one over the size of source box)
  void M_to_T(target_t *first, target_t *last, double scale) const;

  // TODO
  /// Apply a local expansion to a set of targets
  ///
  /// \param first - the first target point
  /// \param last - one past the last target point
  /// \param scale - scaling factor (one over the size of the target box)
  void L_to_T(target_t *first, target_t *last, double scale) const;

  /// Compute the direct interaction between sources and targets
  ///
  /// \param s_first - the first source point
  /// \param s_last - one past the last source point
  /// \param t_first - the first target point
  /// \param t_last - one past the last target point
  void S_to_T(source_t *s_first, source_t *s_last,
              target_t *t_first, target_t *t_last) const;

  // TODO - make sure these are all; finish writeup
  //
  // I think the plan here is to not have to worry about views for the
  // operations. The various operations get their indices, and if the needed
  // views are not available, it should just have some error. Probably an
  // assertion failure to start.
  void M_to_I(Index s_index) const;
  void I_to_I(Index s_index, Index t_index) const;
  void I_to_L(Index t_index) const;

  /// Add an expansion to this expansion
  ///
  /// Given another expansion (assumed to be of the same type), add the
  /// expansions together.
  ///
  void add_expansion(const expansion_t *temp1);
};

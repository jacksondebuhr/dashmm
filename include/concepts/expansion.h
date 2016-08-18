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


/// First, some background information:
///
/// The concept of Expansion in DASHMM is similar to but distint from the
/// concept of a mathematical expansion of a given kernel function. The
/// DASHMM Expansion object is better thought of as a collection of
/// mathematical expansion. This allows for advances methods to be applied
/// using DASHMM, such as the merge-and-shift technique, which results in
/// sets of expansions for a given node of the tree.
///
/// The first distinction between mathematical expansions and the concept
/// in DASHMM is that the Expansion objects will serve multiple sorts of
/// expansions during a typical calculation. This notion is represented as
/// the ExpansionRole, a required argument for constructing most Expansion
/// objects. The mathematical nature of the expansion is often different for
/// nodes in the Source tree and nodes in the Target tree. Further, there
/// are intermediate expansions that might be employed in advanced versions
/// of multipole methods. These intermediate expansions may also take different
/// forms if they are associated with a Source node or a Target node.
///
/// Second, to allow for the possibility that one node of the DAG, that is,
/// the part of the compuation represented by an Expansion, might need multiple
/// mathematical expanions. One use for this is to enable certain advanced
/// multipole methods, but one might also use the same evaluation to perform
/// multipole method computations on multiple kernels at the same time.
///
/// Thus, Expansion implementations are required to present a set of views
/// of the contained data. Each view is one mathematical expansion, each of
/// which is logically related in some way, be that by supplying different
/// required expansions for an advanced method, or by representing other
/// forces involved in a given problem.
///
/// Not all operations for an Expansion will need every view to be performed,
/// and not all views will be created by any given operation. The exact
/// details of which operations need which views are not the concern of
/// DASHMM, but rather are the concern of the implementor of the specific
/// Expansion class.


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
  /// requested. The listed scale is when the expansion type is not
  /// scale-invariant (e.g. Yukawa) Further, the provided role will indicate
  /// how the Expansion will be used. There is no requirement to use any of
  /// these arguments, but the constructor must accept them.
  ///
  /// The scale value that will be provided to the expansion will be provided
  /// by the compute_scale static member of this class.
  Expansion(Point center, double scale, ExpansionRole role);

  /// The second creates the expansion from previously existing data specified
  /// with a ViewSet object.
  ///
  /// Further, this constructor needs to be able to operate in a 'shallow'
  /// mode, where views contains no views. This allows for situations where
  /// Expansion is needed, but the specific data for the expansion is not.
  /// The exemplar of this use is to perform an S->T operation. The expansion
  /// will get a value for n_digits from the provided ViewSet.
  Expansion(const ViewSet &views);

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
  /// release() will just set this object's member to nullptr.
  void release();

  /// Returns if the indicated views are valid
  ///
  /// An expansion is valid if it has data associated with it. After calling
  /// release(), an expansion will be invalidated.
  ///
  /// If view is empty, this will check all views.
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

  /// Get the role of this expansion
  ExpansionRole role() const;

  /// The point around which the expansion is defined.
  Point center() const;

  /// The number of terms in the expansion for the specified view.
  ///
  /// Do not confuse this with the size of the data in bytes for the given
  /// view. That may be retried using get_views().
  size_t view_size(int view) const;

  /// Get a term of the expansion.
  ///
  /// The input should be in the range [0, size()). The particular
  /// implementation will decide if this is range checked. To be safe, assume
  /// that the input is not range checked.
  ///
  /// \param view - the view in question
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
  ///
  /// \returns - The resulting multipole expansion
  std::unique_ptr<expansion_t> S_to_M(Point center, source_t *first,
                                      source_t *last) const;

  /// Create a local expansion for a given set of source points
  ///
  /// This uses the given sources to create a local expansion centered
  /// at \p center. The sources are provided as pointers to the \p first and
  /// one past the \p last source record.
  ///
  /// \param center - the point around which to form the expansion
  /// \param first - address of the first source
  /// \param last - address of one past the last source
  ///
  /// \returns - The resulting local expansion
  std::unique_ptr<expansion_t> S_to_L(Point center, source_t *first, 
                                      source_t *last) const; 

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

  /// Apply a multipole expansion to a set of targets
  ///
  /// \param first - the first target point
  /// \param last - one past the last target point
  void M_to_T(target_t *first, target_t *last) const;

  /// Apply a local expansion to a set of targets
  ///
  /// \param first - the first target point
  /// \param last - one past the last target point
  void L_to_T(target_t *first, target_t *last) const;

  /// Compute the direct interaction between sources and targets
  ///
  /// \param s_first - the first source point
  /// \param s_last - one past the last source point
  /// \param t_first - the first target point
  /// \param t_last - one past the last target point
  void S_to_T(source_t *s_first, source_t *s_last,
              target_t *t_first, target_t *t_last) const;

  /// Compute an intermediate expansion from a multipole expansion.
  ///
  /// \param s_index - the index of the tree node for which this expansion
  ///                  applies
  ///
  /// \returns - the resulting intermediate expansion.
  std::unique_ptr<expansion_t> M_to_I(Index s_index) const;

  /// Compute an intermediate expansion from an intermediate expansion.
  ///
  /// \param s_index - the index of the tree node for which this expansion
  ///                  applies
  /// \param s_size - the size of the tree node for which this expansion
  ///                 applies
  /// \param t_index - the index of the tree node for which the resulting
  ///                  expansion should apply.
  ///
  /// \returns - the resulting intermediate expansion.
  std::unique_ptr<expansion_t> I_to_I(Index s_index, double s_size,
                                      Index t_index) const;

  /// Compute a local expansion from an intermediate expansion.
  ///
  /// \param t_index - the index of the tree node for which this expansion
  ///                  applies
  /// \param t_size - the size of the tree node for which this expansion
  ///                 applies.
  ///
  /// \returns - the resulting local expansion.
  std::unique_ptr<expansion_t> I_to_L(Index t_index, double t_size) const;

  /// Add an expansion to this expansion
  ///
  /// Given another expansion (assumed to be of the same type), add the
  /// expansions together.
  ///
  void add_expansion(const expansion_t *temp1);


  /// Update a kernel table
  ///
  /// This should generate or update a kernel table for the expansion type on
  /// which this was called. These tables have implementations that are
  /// entirely up to the implementer. Typically, these would be the one time
  /// computations that are needed for the expansion in question. The input
  /// give the various factors that might determine the values in the table.
  /// It is likely best if the implementation saves these values, but it is
  /// not required.
  ///
  /// \param n_digits - the number of digits of accuracy required by the
  ///                   expansion
  /// \param domain_size - the size of the top level domain
  /// \param kernel_params - the kernel parameters that specify the interaction
  static void update_table(int n_digits, double domain_size,
                           const std::vector<double> &kernel_params);

  /// Delete a kernel table
  ///
  /// This will delete all tables, if they should exist, associated with the
  /// type of expansion.
  static void delete_table();

  /// Compute the scale to pass into expansion constructors
  ///
  /// This will only be called after the table exists, so the implementation
  /// should rely on the existence of the table. In particular, any parameters
  /// to the kernel should likely be stored in the table, and so if there is
  /// a need to use a kernel parameter that can be assumed to be stored in the
  /// table. Of course, the implementation of the table is responsible for
  /// saving the kernel parameters. Further, if the scale relies on the domain
  /// size in its computation, that too should be saved in the table.
  ///
  /// \param index - the index of the node containing this expansion
  ///
  /// \returns - the appropriate scale factor to be used in an expansion's
  ///            constructor
  static double compute_scale(Index index);
};

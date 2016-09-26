# Basic Guide to DASHMM #

In this chapter, the basic interface to DASHMM will be covered. Generally
speaking, the basic user interface to DASHMM is anything needed to employ the
provided methods and kernels in applications. For instructions on how to define
and use new methods and kernels, please see the next chapter on the advance
user interface.

One of the goals of DASHMM is to make it easy to use multipole methods, and so
the basic interface targets ease-of-use. Additionally, the library aims to make
it easy to perform parallel computations using the multipole method. However,
the advanced dynamic adaptive techniques that DASHMM employs are not simple
to use directly, so the library is constructed in a way that allows the user
to get the benefits of parallel execution without having to write the parallel
code themselved. That being said, it can be useful to have a sense of what
underlies the conceptual framework of DASHMM. And so, the following section will
cover these concepts at a level that might be relevant for basic use of the
library. More details can be found in the next chapter.

## DASHMM Concepts ##

This section covers both the conceptual framework of the multipole methods
supported by DASHMM and the parallelization of those methods. More explanation
of DASHMM's conceptual framework can be found in the following chapter, or
in the code paper: (Put in citation here once available).

### Multipole Method Abstractions ###

DASHMM is a templated library allowing for the description of a general set of
multipole methods with a small set of template parameters. To use DASHMM, one
must specify four types: the Source data type, the Target data type, the
Expansion type and the Method type. Within certain limits, each of these
can be varied independently, for example, a given expansion might be used with
a number of methods, allowing users to easily experiment and select the method
that most meets their needs.

DASHMM includes a number of built-in expansion and method types. For more, see
below.

#### Source ####

The Source type gives the structure of the source point data. There are few
requirements on the Source type. Typically the minimum requirements are the
position and charge of the source.  The `position` is of type `Point` and the
`charge` is of type `double`. For example, the following is a minimal Source
type that works with every expansion provided with DASHMM.

~~~{.cc}
struct SourceData {
  Point position;
  double charge;
};
~~~

Beyond these required members, anything might be added to a Source type. This
allows the user to associate application specific data to the sources in the
evalution. The only additional requirement is that the resulting type be
trivially copyable.

#### Target ####

The Target type gives the structure of the target point data. There are
similarly few requirements on the Target type. Each Target type will need
a `position` member of type `Point`. Then, each expansion requires a
particular member to store the result. The details on the exact requirements
can be found in the documentation of the individual expansions below.
Typically this is a member `phi` of type `std::complex<double>`, which has
been aliased as `dcomplex_t` in DASHMM. So the following would work for many
of the included DASHMM expansions:

~~~{.cc}
struct TargetData {
  Point position;
  std::complex<double> phi;
};
~~~

Beyond the required members, anything might be added to the Target type.
A typical choice is an identifier of some kind because DASHMM evaluations will
sort the input data to suit the parallel computation, and points can then be
identified after the computation. Like with the Source type, the resulting
Target type must be trivially copyable.

#### Expansion ####

The particular potential or interaction that is being computed with the
multipole moment is called the kernel. For example, the Laplace kernel is the
traditional potential from electrostatics or Newtonian gravitation. In
DASHMM, kernels are not represented directly. Instead, Expansions are
created that implement the needed operations for the particular kernel. The
distinction is that there are often multiple ways to expand a given kernel to
be used in a multipole method computation. Each Expansion represents a way
(or a closely related set of ways) that the potential is expanded into an
approximation.

The Expansion type is a template type over two paramters, the source and
target types being employed. It is the Expansion type that places restrictions
on the Source and Target types; the Expansion requires certain data to
compute from (taken from the Sources), and it will produce certain data (into
the Targets). This allows the Expansion to operate on any data that meets its
requirements, meaning that a very general set of source and target data types
can be supported.

For details on creating user-defined Expansions, please see the advanced use
guide in the next chapter. This will also include the set of operations that
an Expansion must support.

#### Method ####

The final major abstraction in DASHMM is the Method. This type specifies how
the Expansion is used on the provided Source data to compute the interaction
at the locations specified in the Target data. The Method is responsible for
connecting the various Expansions representing the hierarchically subdivided
set of Source and Target locations with the appropriate operations provided
by the Expansion. It is the Method that allows one to perform both a Barnes-Hut
computation as well as the Fast Multipole Method.

The Method type is a template over three types: the Source, the Target and the
Expansion. For the Method, it is not important exactly how the various
operations are implemented, but only that they are implemented. The details of
the Expansion, and thus the central details of the particular interaction
being studied, are hidden and unimportant to the method.

Despite this generality, it is possible to create a Method that only works for
certain expansions. Indeed, included in DASHMM is the `LaplaceCOM` expansion,
which does not provide implementations for all of the operations needed by
the `FMM` or `FMM97` methods. This is intentional, as the style of expansion
in `LaplaceCOM` is not terribly well suited to the Fast Multipole Method.
Nevertheless, great flexibility and generality is possible with DASHMM.

For details on creating user-defined Methods, please see the advanced use guide
in the next chapter. This will also include a description of the set of
routines that make up the method concept.

### Parallelization Abstractions ###

DASHMM uses the advanced runtime system HPX-5 for its parallelization. HPX-5
provides a number of features that allow DASHMM to naturally express the
parallelism and data dependence of the computation directly in programming
constructs. However, the flexibility of HPX-5 comes with a significant amount
of effort to learn the system. One major goal of DASHMM is to get the benefits
of the dynamic adaptive techniques enabled by HPX-5 without the end-user having
to write directly to HPX-5 constructs, and DASHMM is successful in this regard.
It is, nevertheless, useful to have some notion of a few concepts from HPX-5
for the basic use of DASHMM. For more details on HPX-5, please visit
https://hpx.crest.iu.edu/. For more details on the use of HPX-5 in DASHMM,
please see the advanced use guide in the next chapter, or the paper FILL IN
DETAILS OF PAPER.

#### PGAS ####

HPX-5 provides a partitioned global address space (PGAS) that DASHMM uses for
the data during the computation.

That the address space is partitioned means that though the addresses are
unified into a global space, each byte of the global address space is served
by the physical memory on one particular locality of the system. A locality
is similar to the concept of a rank in an MPI program. That the address space
is global means that there is a single virtual address space allowing any
locality to refer to data, even if that data is not stored in the same physical
memory of the referring locality.

Practically speaking for users of DASHMM, the important part of the global
address space provided by HPX-5, as used by DASHMM, is that it is partitioned.
Each locality will have direct access to a portion of the data, and more
indirect access to all of it.

#### Execution Model ####

The execution model for a DASHMM program is one very similar to the SPMD model.
The program written to use DASHMM will be run on every locality in the allocated
resources. And at certain points, the exection will be handed off to DASHMM and
ultimately HPX-5 by making DASHMM library calls. Inside DASHMM, the execution is
very dynamic and is formed of a large number of small, interdependent tasks.
This complication, however, is hidden from the user. Instead, DASHMM presents
an interface where each locality participates in collective calls, with each
presenting possibly different data to the DASHMM library. So in many ways,
using basic DASHMM will be very similar to using MPI.

## Basic Types ##

DASHMM defined a number of basic types that are used throughout the system,
and which might be used by users of the library.

### `ReturnCode` ###

DASHMM calls will return values of the type `ReturnCode` where it is reasonable
to do so. The possible values are: `kSuccess`, `kRuntimeError`, `kIncompatible`,
`kAllocationError`, `kInitError`, `kFiniError` and `kDomainError`. See specific
library calls for cases where each might be returned.

### `dcomplex_t` ###

Many kernels return potential values that are complex numbers. DASHMM provides
`dcomplex_t` an alias to `std::complex<double>`.

### `Point` ###

The `Point` class is used to represent locations in three dimensional space.
`Point` is expected for giving the locations of sources and targets. It has the
following members:

 - `Point::Point(double x = 0, double y = 0, double z = 0)` : Construct a
 point from a given set of coordinates. This can also be used to default
 construct a point.
 - `Point::Point(double *arr)` : Construct a point from a C-style array. It is
 important that `arr` should contain at least three members.
 - `Point::Point(const Point &pt)` : Copy construct a point.
 - `Point Point::scale(double c) const` : Return a point whose coordinates have
 all been scaled by the factor `c`.
 - `double Point::operator[](size_t i) const` : Indexing access to the
 coordinates of the point. `i` must be in the range `[0,2]`.
 - `double Point::x() const` : Return the x coordinate of the point.
 - `double Point::y() const` : Return the y coordinate of the point.
 - `double Point::z() const` : Return the z coordinate of the point.
 - `double Point::norm() const` : Return the 2-norm of the point.
 - `void Point::lower_bound(const Point &other)` : This computes the lowest
 coordinate in each direction of this point and `other` and sets this point's
 coordinate to that value.
 - `void Point::upper_bound(const Point &other)` : This computes the highest
 coordinate in each direction of this point and `other` and sets this point's
 coordinates to that value.

Additionally, there are a few non-member operations that are defined:

 - `double point_dot(const Point &left, const Point &right)` : Treat the
 points as if they are vectors and take their dot product.
 - `Point point_add(const Point &left, const Point &right)` : Perform a
 component-wise addition of `left` and `right` and return a point with
 the result.
 - `Point point_sub(const Point &left, const Point &right)` : Perform a
 component-wise subtraction of `right` from `left` and return a point with
 the result.

## Initializing DASHMM ##

DASHMM must be initialized and finalized to be used. There are some DASHMM
operations that must only occur before initialization, and some that can only
occur after initialization. All DASHMM operations must occur before the
library is finalized.

### `ReturnCode init(int *argc, char ***argv)` ###

Initialize the runtime system supporting DASHMM and allocate any resources
needed by DASHMM. The addresses of the command line arguments must be
provided as the behavior of HPX-5 can be controlled by these arguments. Any
arguments dealing with HPX-5 directly will be removed and `argc` and `argv`
will be updated accordingly.

`init()` returns `kSuccess` if the system is successfully started, and
`kRuntimeError` otherwise. If `init()` returns `kRuntimeError` all subsequent
calls to DASHMM will have undefined behavior.

All other DASHMM library calls must occur after `init()`. However, some DASHMM
related objects must be constructed before the call to `init()`. See below for
details.

This is a collective call; all localities must participate.

### `ReturnCode finalize()` ###

This will free any DASHMM specific resources and shut down the runtime system.
No other calls to DASHMM must occur after the call to `finalize()`. This is
a collective call; all localities must participate.

## SPMD Utilities ##

DASHMM provides a small number of traditional SPMD utilities to make certain
things simpler.

### `int get_my_rank()` ###

This returns the rank of the calling locality.

### `int get_num_ranks()` ###

This returns the number of ranks available.

### `void broadcast(T *value)` ###

This performs a broadcast of the given value at rank 0 to all other ranks.
This is a template over the type `T`. This is a collective operation. Each rank
must provide the address of a type `T` object. For rank 0, this will provide
the address of the value to share; for all other ranks this provides the address
into which the value broadcast from rank 0 will be stored.

## Evaluation ##

The central object in any DASHMM evaluation is the `Evaluator` object. This
object not only manages the registration of certain actions with the runtime
system, but also provides the interface to performing the multipole method
evaluation.

The `Evaluator` object is a template over four types: the source type, the
target type, the expansion type and the method type. The `Evaluator` for a
given set of types must be declared before the call to `init()`. For example:

~~~{.cc}
dashmm::Evaluator<Source, Target, dashmm::Laplace, dashmm::FMM97> eval{};
~~~

would be an `Evaluator` for the Laplace kernel using the advanced FMM method
(see below for details) for two user-defined types `Source` and `Target`
implementing the data that the user requires of the source and target points.

Specifying the full type of the evaluator will cause the template to expand
out all of the needed actions to actually implement the evaluation using
HPX-5, and will also register those actions with HPX-5.

### `ReturnCode Evaluator::evaluate(...)` ###

Perform a multipole method evaluation. The arguments to this method, in order,
are as follows:

 - `const Array<Source> &sources` : the Array containing the source data.
 - `const Array<Target> &targets` : the Array containing the target data.
 - `int refinement_limit` : the refinement limit of the tree. The sources and
 targets will be placed into a hierarchical partitioning of space. This
 partitioning will end when there are fewer sources or targets than the supplied
 refinement limit.
 - `const Method<Source, Target, Expansion<Source, Target>> &method` : the
 method to use for the evaluation. A few methods require parameters at
 construction, so this is passed in to provide those paramters to the
 evaluation.
 - `int n_digits` : the accuracy parameter for the `Expansion` in use. This is
 often the number of digits of accuracy required.
 - `const std::vector<double> &kernelparams` : the parameters for the kernel
 evaluation. These are those quantities that are constant for each use of the
 particular expansion. For example, the strength of the electrostatic force
 when using the Laplace kernel is one example of a kernel parameters. See the
 individual expansions for details about what needs to be provided.

Note that during evaluation, the records in the source and target arrays may be
sorted. So a separate identifier should be added to the source and target types
if the identity of the sources or targets needs to be tracked. Other than the
sorting, the only change to the data in the target array will be the output
potential or other field value (as specified by the chosen expansion). The
only change to the source array beyond the sorting will only occur in the case
that the source and target arrays are the same, in which case the previous
comment about the targets also applies to the source.

This is a collective call; all localities must partiticipate.

The possible return values are `kSuccess` when there is no problem, and
`kRuntimeError` when there is a problem with the execution.

## DASHMM Array ##

## Built-in Methods ##

## Built-in Expansions ##
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
code themselves. That being said, it can be useful to have a sense of what
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
evaluation. The only additional requirement is that the resulting type be
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

The Expansion type is a template type over two parameters, the source and
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
resources. And at certain points, the execution will be handed off to DASHMM and
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
 construction, so this is passed in to provide those parameters to the
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

This is a collective call; all localities must participate.

The possible return values are `kSuccess` when there is no problem, and
`kRuntimeError` when there is a problem with the execution.

## DASHMM Array ##

DASHMM provides an array construct that represents a collection of records in
the global address space provided by HPX-5. The `Array` object is a template
type over the record. The only requirement on the record type is that it is
trivially copyable. When an `Array` is created, as in

~~~{.cc}
Array<T> source_data{};
~~~

no global memory is yet allocated for the array. To allocate the global memory
that will serve the array, one must use `allocate()` (see below for details).
Array objects can be thought of as a traditional array, but which is broken
into one part for each locality in the system. The parts, or segments of the
array, could have different lengths. Some segments can even have no records,
meaning that a given locality does not have a share of the data.

`Array` objects are intended to be employed in user code, and so most members
of the interface are SPMD in nature, and are collective operations.

### `Array<T>::Array()` ###

`Array` objects have a single constructor, which does not allocate any global
memory. It takes a single optional argument that is not needed for basic use
of DASHMM. When default constructed, an `Array` will be invalid (see `valid()`
below.)

This is not a collective operation. To use an array, each locality must create
an `Array` object. The individual `Array` objects will be bound into a unified
object referring to the global address space using `allocate()`.

### `bool Array<T>::valid() const` ###

This returns if the array is valid, that is, it refers to some global memory.
When initially constructed, an `Array` will be invalid, and can only be made
valid by allocating global memory for the records (see `allocate()`).

This method is not collective.

### `size_t Array<T>::count() const` ###

This returns the number of records in the segment of the array on the calling
locality. This is a collective call, and it is an error to call this method
on an invalid array. Each calling rank will receive a different result from
this method.

### `size_t Array<T>::length() const` ###

This returns the total length of the entire array. The result of `length()` is
equal to the sum of the results of `count()` from each rank. This is a
collective call, and it is an error to call this method on an invalid array.

### `ReturnCode Array<T>::allocate(size_t record_count)` ###

This allocated the global memory to serve an array with the given counts on
each locality. After this call the array will be valid unless there is an
error, which will be indicated by the return code. Possible return values are:
`kSuccess` if the array is successfully allocated; `kDomainError` is the object
already has an allocation; `kAllocationError` is the global memory cannot be
allocated; or `kRuntimeError` if there is some error in the runtime.

This is a collective call. Each locality can provided a different number of
records to allocate via the `record_count` parameter. The resulting global
allocation will match the input of `allocate()`. A locality can provide `0`
as an argument, so long as at least one rank asks for a non-zero number of
records. For instance, one locality might allocate all of the records for
an array, or each locality might own a portion of the overall records.

### `ReturnCode Array<T>::destroy()` ###

This will destroy the global memory allocated for this array object. It is
an error to call this on an invalid array. This is a collective call. The
result of the method will be either `kRuntimeError` if there is an error
in the runtime, or `kSuccess` otherwise.

### `ReturnCode Array<T>::get(size_t first, size_t last, T *out)` ###

This method gets data from an array object, and places the requested
records into the buffer provided by `out`. The range of records that is
retrieved is specified by `first` (inclusive) and `last` (exclusive).

This is a collective call. Each locality will provide a different range of
records, and a different local buffer into which the retrieved data will be
placed. `first` and `last` are given in terms of the segment of the array on
this locality. To discover the number of records on the calling locality, use
`count()`.

Note that this is a copy of the global data; changes to the retrieved values
will not be reflected in the array object.

This method will return one of the following: `kRuntimeError` if there is an
error with the runtime; `kDomainError` if the provided index range is
inconsistent with the array object; or `kSuccess` otherwise.

### `ReturnCode Array<T>::put(size_t first, size_t last, T *in)` ###

This method puts data into an array object, copying the specified
records from the buffer provided by `in`. The range of records that is
copied is specified by `first` (inclusive) and `last` (exclusive).

This is a collective call. Each locality will provide a different range of
records, and a different local buffer from which the retrieved data will be
placed. `first` and `last` are given in terms of the segment of the array on
this locality. To discover the capacity of the array on the calling locality,
use `count()`.

Note that this places a copy of the records into the global address space;
changes in the local data will not be reflected in the array object.

This method will return one of the following: `kRuntimeError` if there is an
error with the runtime; `kDomainError` if the provided index range is
inconsistent with the array object; or `kSuccess` otherwise.

### `T *Array<T>::collect()` ###

This method collects all of the records in the array and returns a new local
allocation containing the records at locality 0. This is largely speaking a
convenience feature. Note that this will allocate a copy of the entirety of the
array on one locality.

This is a collective call. The returned pointer will only be valid on locality
zero. All other ranks will receive `nullptr`.

Note that this provides a copy of the data; changes to the returned data will
not be reflected in the array object.

## Array Map Actions ##

To avoid the round trip from and to the global address space via `Array`'s
`get()` and `put()` methods, one can make use of the `ArrayMapAction` type.
This type specifies an action to be performed on the records of an array.
Then, together with a method of the `Array` object, this allows for
some computation to occur on the records of an array.

This class is a template requiring two parameters: the type of records for
the array to which the action will be applied (hereafter `T`), and a type
specifying an environment to provide to the action (hereafter `E`).

The action is specified during construction of an object of type
`ArrayMapAction`. This is provided as a function pointer, for a function with
a particular signature. So that DASHMM can use the action on every locality,
like the `Evaluator` object, any `ArrayMapAction` objects must be defined
before `init()` is called.

The function implementing the action must take three arguments, the first is
a `T *` giving the data on which to act, the second is a `const size_t` gives
the number of records on which to act, and the third is `const E *`, where `E`
is the environment type of the `ArrayMapAction`. For example:

~~~{.cc}
void update_position(T *data, const size_t count, const E *env) {
  for (size_t i = 0; i < count; ++i) {
    data[i].position += data[i].velocity * env->delta_t;
  }
}
~~~

is an action that might perform a position update for a time-stepping code.
It is important to note that the action implementation can place requirements
on the array record type `T`. In the previous example, type `T` needs to have
a member called `velocity`. The first argument to the action will be a
correctly offset pointer into a contiguous chunk of records. The action should
only assume that `data[0]` through `data[count - 1]` are available for use.

### `ArrayMapAction<T, E>::map_function_t` ###

This class aliases the type of function that can implement the action. It is:

~~~{.cc}
void (*)(T *, const size_t, const E *)
~~~

### `ArrayMapAction<T, E>::ArrayMapAction(map_function_t f)` ###

This constructs an array map action object. The provided function pointer
will give the action that is performed on the array. To use an action,
the associated `ArrayMapAction` must be defined before `init()` is called.
This object will register the needed actions with the runtime system.

For any given action `f` only a single `ArrayMapAction` should be defined.

### `ReturnCode Array<T>::map(ArrayMapAction<T, E> &act, const E *env)` ###

Once an `ArrayMapAction` is defined, it can be used on a specific array by
calling that array's `map()` method. This will cause the action represented
by `act` to be applied on all entrieds of this array. The action ultimately
works on segments of the array. The environment, `env` is provided unmodified
to each segment.

This is a collective call. It is an error to `map()` on an invalid array.

## Built-in Methods ##

DASHMM includes a number of built-in methods that are ready to use for
problems: the Barnes-Hut method, two forms of the Fast Multipole Method and
a Direct summation method. These will be covered in detail below. To
successfully compile, all operations in the Expansion concept need to be
implemented for a given expansion. However, some methods only require a
subset of the full complement of operations. These will be covered below,
in order of increasing complexity.

### `Direct` ###

The `Direct` method is primarily intended for use as a comparison for computing
the exact result for a given set of sources. There is some parallelism
implemented for this method, so the cost can be reduced somewhat. However, for
realistic problem sizes, this method should be avoided.

The only operation that is used by the `Direct` method is the `S->T` operation,
so any Expansion implementing that operator can be used with the `Direct`
method.

Please see the demo program `demo/basic` included with DASHMM for an example
use case of `Direct` and note that it is only applied to a very small subset of
the target locations.

Construction of `Direct` is simple, as it can only be default constructed.

### `BH` ###

The `BH` method implements the Barnes-Hut algorithm in the framework of
DASHMM. The implemented method uses the simple multipole acceptance criterion
parameterized by a single angle. If a multipole expansion is centered a
distance `R` away from a given point and the multipole expansion is associated
with a region of size `D` then the multipole expansion is used only if
`D/R < theta_C` for some critical angle `theta_C`.

In practice, since the
hierarchical partioning of the source and target locations in DASHMM produces
leaves that can contain multiple particles, the multipole acceptance criterion
is evaluated pessimistically: the radius is computed from the nearest possible
location in a given target tree node to the multipole in question. This will
tend to produce results with a slightly smaller error, with a slightly larger
execution time.

When creating a `BH` object, the opening angle `theta_C` is specified:

~~~{.cc}
BH bh_method(0.6);
~~~

Generally the opening angle should be less than one. As the angle approaches
zero, the method approaches the direct summation method. The critical angle
for a particular `BH` object can be accessed with `double BH::theta()`,
which returns the angle specified at creation time. A default constructed
`BH` uses a critical angle of `0`.

In addition to the direct contribution, `S->T`, to use `BH` an expansion must
also implement fully the following operations: `S->M`, `M->M` and `M->T`.

### `FMM` ###

The `FMM` method implements the Fast Multipole Method in its original form.
To create an `FMM` method, no arguments are needed as the decisions about
which expansions to use in which situations are all made based on fixed
geometric considerations.

The following operations must have a full implementation in an expansion to be
used with `FMM`: `S->T`, `S->M`, `S->L`, `M->M`, `M->L`, `M->T`, `L->L`,
and `L->T`.

### `FMM97` ###

The `FMM97` method implements the Fast Multipole Method in the form that
uses exponential expansions and the merge-and-shift technique. No parameters
are needed to construct an `FMM97` method.

The following operations must have a full implementation in an expansion to be
used with `FMM97`: `S->T`, `S->M`, `S->L`, `M->M`, `M->L`, `M->T`, `L->L`,
`L->T`, `M->I`, `I->I` and `I->L`.

## Built-in Expansions ##

DASHMM includes a number of built-in expansion that are ready to use for
applications. Each expansion will have a different set of implemented
operations which will restrict their use for certain methods. Further, each
expansion will place requirements on the source and target types used during
a DASHMM evaluation. These details will be covered for each expansion
below.

### `Laplace` ###

The `Laplace` expansion is a spherical harmonic expansion of the Laplace
potential. Unlike `LaplaceCOM` and `LaplaceCOMAcc` below, this expansion is
designed to handle sources with both signs of charge without losing accuracy.
This potential is scale-invariant, so there are no kernel parameters that are
needed in the call to `evaluate()`. However, calls to evaluate must supply
an accuracy parameter giving the number of digits of accuracy that are
requested.

Though this expansion is in principle compatible with every method included
with DASHMM, it is designed for the `FMM` and `FMM97` methods.

This expansion imposes the following restrictions on the source type: a
member of type `Point` with the name `position` must be provided; a member of
type `double` with the name `charge` must be provided.

This expansion imposes the following restrictions on the target type: a
member of type `Point` with the name `position` must be provided; a member of
type `dcomplex_t` with the name `potential` must be provided.

### `Yukawa` ###

The `Yukawa` expansion is a spherical harmonic expansion of the Yukawa
potential. This potential is scaling variant, so a single kernel parameter
must be provided to `evaluate()`. Calls to evaluate must also supply an
accuracy parameter that gives the number of digits of accuracy required.

Though this expansion is in principle compatible with every method included
with DASHMM, it is designed for the `FMM` and `FMM97` methods.

This expansion imposes the following restrictions on the source type: a
member of type `Point` with the name `position` must be provided; a member of
type `double` with the name `charge` must be provided.

This expansion imposes the following restrictions on the target type: a
member of type `Point` with the name `position` must be provided; a member of
type `dcomplex_t` with the name `potential` must be provided.

### `LaplaceCOM` ###

The `LaplaceCOM` expansion is a center of mass expansion of the Laplace
potential. This form of the expansion extends to the quadrupole term, and
because of the choice of center has an identically zero dipole term. This
expansion can be used to compute the potential. See `LaplaceCOMAcc` for an
equivalent expansion that computes the acceleration. This potential is
scale-invariant so not kernel parameters are needed in the call to `evaluate()`.
Further, the number of terms in the expansion is fixed, so the accuracy
parameter to `evaluate()` is ignored.

This expansion is only compatible with the `BH` or `Direct` methods; it does
not implement all the needed operations for use with `FMM` or `FMM97`.

This expansion imposes the following restrictions on the source type: a
member of type `Point` with the name `position` must be provided; a member of
type `double` with the name `charge` must be provided.

This expansion imposes the following restrictions on the target type: a
member of type `Point` with the name `position` must be provided; a member of
type `dcomplex_t` with the name `potential` must be provided.

### `LaplaceCOMAcc` ###

The `LaplaceCOM` expansion is a center of mass expansion of the Laplace
potential. This form of the expansion extends to the quadrupole term, and
because of the choice of center has an identically zero dipole term. This
expansion can be used to compute the acceleration. See `LaplaceCOM` for an
equivalent expansion that computes the potential. This potential is
scale-invariant so not kernel parameters are needed in the call to `evaluate()`.
Further, the number of terms in the expansion is fixed, so the accuracy
parameter to `evaluate()` is ignored.

This expansion is only compatible with the `BH` or `Direct` methods; it does
not implement all the needed operations for use with `FMM` or `FMM97`.

This expansion imposes the following restrictions on the source type: a
member of type `Point` with the name `position` must be provided; a member of
type `double` with the name `charge` must be provided.

This expansion imposes the following restrictions on the target type: a
member of type `Point` with the name `position` must be provided; a member of
type `double [3]` with the name `acceleration` must be provided.

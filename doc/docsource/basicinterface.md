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

- parallel through HPX-5, but hidden
- PGAS
- overall execution model SPMD with HPX, or single rank diffusive

## Basic Types ##

## Initializing DASHMM ##

## Evaluation ##

## DASHMM Array ##

## Built-in Methods ##

## Built-in Expansions ##
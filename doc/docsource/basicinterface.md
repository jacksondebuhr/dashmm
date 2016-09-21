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

- expansion, method, source target
- parallel through HPX-5, but hidden
- PGAS
- overall execution model SPMD with HPX, or single rank diffusive

## Basic Types ##

## Initializing DASHMM ##

## Evaluation ##

## DASHMM Array ##

## Built-in Methods ##

## Built-in Expansions ##
# Introduction to DASHMM #

The Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM) is a
C++ library providing a general framework for computations using multipole
methods. In addition to the flexibility to handle user-specified methods and
expansion, DASHMM includes built-in methods and expansions, including the
Barnes-Hut (BH) and Fast Multipole Method (FMM), and expansions
implementing the Laplace and Yukawa kernels.

DASHMM is designed to make its adoption and use as easy as possible, and so the
interface to DASHMM provides both an easy-to-use basic interface and a more
advanced interface. The basic interface allows a user to get up and running
with multipole methods as quickly as possible. However, more advanced use-cases,
especially for users that are implementing their own methods and expansions,
will need to explore the advanced interface. Finally, for the truly interested,
the entire library reference documentation can be found at the end of this
document.

DASHMM is built using the advanced runtime system, HPX-5, but the basic, and
much of the advanced interface does not require any knowledge of how to use
HPX-5 directly. Instead, DASHMM insulates the users from the specific details
of HPX-5, allowing expression of the multipole method application in higher
level concepts than the specifics of threads and execution control structures.
Nevertheless, it can be helpful to have a sense of the conceptual underpinnings
of HPX-5, and how those relate to DASHMM. Information about this can be found
in this guide in both the basic and advanced interface. Even more information
can be found at the HPX-5 website (https://hpx.crest.iu.edu/).

TODO: add a note about the first code paper.

This document covers version 1.0.0 of DASHMM. For the latest news and updates
on DASHMM, please visit the DASHMM website:
https://www.crest.iu.edu/projects/dashmm/.

The DASHMM project has adopted Semantic versioning. For a description,
please see http://semver.org/.

In the following, snippets of code or the names of code constructs will be
set in a fixed width font. For example, `main()`. Unless otherwise indicated,
every construct presented in this guide is a member of the `dashmm`
namespace.


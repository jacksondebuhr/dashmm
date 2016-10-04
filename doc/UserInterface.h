// This is a helper file to catalog all of the elements of the intended user
// inferface to DASHMM. This will aid in the creation of the user guides
// and the reference documentation. Further, this will allow for easier
// rearrangement of the source, as this will collect everything that a user
// might be expected to use.

// what is the difference between the basic and advanced user guide? The
// advanced has the more interesting parts of the user interface.
// This means there likely should be a third thing, a complete reference
// document. That latter one should likely be built automatically from the
// source. This might be helpful for other things as well. Using groups and
// so on. We can likely even put the rest of the documentation into the
// Doxygen markdown format, and then just have it automatically compile it into
// the various formats that we will need.
//
// Another difference is that some things become relevant if one is creating
// a new expansion or method. The Advanced Guide should cover the needed
// information in those cases as well. e.g. ViewSet


// arraymapaction.h
ArrayMapAction::ArrayMapAction()
//--- did not cover the final template argument


/*


Basic User Guide Outline


I. Intro
II. Installation
 - prerequisites
 - installation
 - using dashmm
 - test codes
III. Overview of Important Concepts
 - expansion, method, source target
 - parallelism though HPX-5, but hidden away
 - partitioned global address space
 - overall execution model SPMD with some HPX, or single rank diffusive
IV. Basic Types
V. init and fini
IX. Evaluation
VI. Array and related
VII. Available Expansions
VIII. Available Methods




*/
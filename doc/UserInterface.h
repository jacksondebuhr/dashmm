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

// array.h
Array::Array()
bool Array::valid()
size_t Array::count()
size_t Array::length()
ReturnCode Array::allocate()
ReturnCode Array::destroy()
ReturnCode Array::get()
ReturnCode Array::put()
ReturnCode Array::map()
T *Array::collect()


// arraymapaction.h
ArrayMapAction::ArrayMapAction()


// evaluator.h
Evaluator::Evaluator()
ReturnCode Evaluator::evaluate()



// types.h -- Are these needed?
ExpansionRole
Operation
// What about DomainGeometry? Index?


// bh_method.h
BH::BH()
double BH::theta()


// direct_method.h
Direct::Direct()


// fmm97_method.h
FMM97::FMM97()


// fmm97distro.h
FMM97Distro::FMM97Distro()


// fmm_method.h
FMM::FMM()


// randomdistro.h
RandomDistro::RandomDistro()


// singlelocdistro.h
SingleLocality::SingleLocality()


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
---
layout: default
title: Tutorial
---

# DASHMM Tutorial

This tutorial covers the use of the basic interface to DASHMM. The basic interface is designed with ease in mind, so that a domain scientist might move rapidly from deciding to use a multipole method to having a working code. One key to this ease is that although DASHMM relies on an advanced runtime system (HPX-5) to provide the parallel computation, the user does not need to know anything about how to use HPX-5 to write a parallel multipole moment utility: The parallelization is handled by DASHMM.

This tutorial will walk through the demonstration utility provided with DASHMM. The code can be found in the `demo/basic/` sub-folder of the DASHMM source code. There is a single source file `testmain.cc` for this tutorial. This tutorial is compatible with the 1.2.0 version of DASHMM.

## Overview of Execution

DASHMM uses the advanced runtime system HPX-5 to manage the parallel computation. HPX-5 has an execution model that is very different from the typical communicating sequential processes (CSP) model that is common in MPI applications. Programs using HPX-5 express their parallelism in terms of a number of lightweight actions, sent across the system using parcels, which are synchronized with a set of lightweight synchronization primitives known as LCOs. However, DASHMM is designed to be easy to use, and so this complication is hidden from the user. From outside, the use of DASHMM will seem to be similar to the typical CSP model: each rank participating in the computation will run the same code, and will at times call out to DASHMM routines collectively. The execution inside the library will be quite different, but this is hidden from the user.

## Getting started

To make DASHMM available for your program, it is sufficient to `#include "dashmm.h"`. This header file will itself include every element of the user interface to DASHMM, including those for advanced use of the library. Every symbol in DASHMM is part of the `dashmm` namespace. In this tutorial, we will explicitly scope all DASHMM constructs to make it clear what is provided by the library, and what must be provided by the user.

We start with `main()`:

```C++
int main(int argc, char **argv) {
  auto err = dashmm::init(&argc, &argv);
  assert(err == dashmm::kSuccess);

  InputArguments inputargs;
  int usage_error = read_arguments(argc, argv, inputargs);

  if (!usage_error) {
    perform_evaluation_test(inputargs);
  }

  err = dashmm::finalize();
  assert(err == dashmm::kSuccess);

  return 0;
}
```

To use DASHMM in a program, one must initialize DASHMM using `dashmm::init()`. This call both initializes the HPX-5 runtime system as well as some internal bookkeeping data used by DASHMM. HPX-5 accepts some command-line arguments that control its behavior, and so we must provide those arguments to DASHMM. HPX-5 will remove the arguments from the command line used by HPX-5, so after the call to `dashmm::init()` the remaining command-line arguments will be ready to be interpreted by the application. For a list of the available HPX-5 related command line arguments, please visit the [HPX-5 website](https://hpx.crest.iu.edu/).

Next this will interpret the incoming arguments, and perform the evaluation test. The first of these is standard, so it will not be covered. The second contains the bulk of the example use of DASHMM, so we shall turn to it after discussing the final call.

`dashmm::finalize()` shuts down the runtime and frees all resources used by DASHMM. There is no reason that this needs to be the last operation in a particular utility. Normal serial work can be performed after the call to `dashmm::finalize()`. However, it is an error to attempt use of any DASHMM routines after `dashmm::finalize()`.

Note also that we are not checking all of the possible return codes from these library calls. Instead, we just `assert()` that everything comes out fine. This will be continued in the rest of the example code. For a list of the possible error codes that can be returned by the library calls, please see the DASHMM Basic User Guide in [Resources]({{ site.baseurl}}/resources.html).

## The Four Main Types

DASHMM uses C++ templates to allow for some flexibility and some easy generalization for the user. For a multipole method evaluation using DASHMM, there are four template parameters that need to be provided. These four parameters give four types that will specialize the general library onto the problem of interest.

The four main types are the Source, Target, Expansion and Method types. There are some light dependencies between these types. The Source and Target will have some restrictions placed on them by the Expansion. The Method will place some restrictions on the Expansion.

For the code in this tutorial, we will have the following Source and Target type:

```C++
struct SourceData {
  dashmm::Point position;
  double charge;
};

struct TargetData {
  dashmm::Point position;
  std::complex<double> phi;
  int index;
};
```

In the previous, the `position`, `charge` and `phi` members are requirements of the expansions that we will be using. This tutorial does not do much other than compute the potential at the targets because of the sources, so these two types do not have too much beyond the required members. We have included an `index` because the sources and targets might be reordered during evaluation, and we will be interested in the original arrangement when we compare the accuracy of the method with the exact result obtained with direct summation.

This is a general feature of the Source and Target types. A user is free to add whatever data they wish to the type, so long as the resulting type is trivially copyable.

## Setting up the Sources and Targets

After some pretty standard argument parsing, which will not be covered here, `main` calls out to `perform_evaluation_test` to perform most of the work of the demo. This routine starts as follows:

```C++
void perform_evaluation_test(InputArguments args) {
  srand(123456);

  dashmm::Array<SourceData> source_handle = prepare_sources(args);
  dashmm::Array<TargetData> target_handle = prepare_targets(args);

  ...
```

DASHMM provides a global address space from which it will read the data on which the evaluations takes place. The `Array<T>` class is a template class over a record type, `T`. This object is a handle to global memory that contains an array of records. These records are shared among the localities in the system, with each locality having some share of the total number of records.  These handles are returned by two utility functions, which are essentially identical. The first of these, `prepare_sources`, is given below.

```C++
dashmm::Array<SourceData> prepare_sources(InputArguments &args) {
  SourceData *sources{nullptr};
  if (args.source_count) {
    sources = new SourceData[args.source_count];
    set_sources(sources, args.source_count, args.source_type, args.method);
  }

  dashmm::Array<SourceData> retval{};
  int err = retval.allocate(args.source_count);
  assert(err == dashmm::kSuccess);
  err = retval.put(0, args.source_count, sources);
  assert(err == dashmm::kSuccess);

  if (args.source_count) {
    delete [] sources;
  }

  return retval;
}
```

In `prepare_sources`, locally allocated `Source` records are filled with randomly generated particle data. The details of the distribution of the points and the associated charges can be controlled with some command line parameters to the demo program. For details on these, please see the `README` included with the demo, or provide `--help` as a command line argument to the demo program.

A DASHMM `Array` is created by calling `allocate()`, providing a number of records. Each rank participating in execution will provide its own number of records to `allocate()`. The result is that the number of records allocated for the `Array` in the global address space will be the sum of all the inputs to `allocate()` across all participating ranks. The distribution of this data will match the inputs: a given rank will own a number of records equal to the argument supplied to `allocate()`. This program is set up so that only rank 0 will have a nonzero value, so all of the records will initially begin on rank 0. This will change during the execution; DASHMM sorts the records to improve the execution time.

Once the global memory is allocated, each rank will put the local memory into the global address space using the `put()` method of `Array`. This copies the local records into the global address space. After this has been accomplished, there is no longer a need for the local data, so it is deleted, and the `Array` is returned.

## Expansions and Methods

Once the data are specified, to perform a multipole method evaluation, one needs to select an Expansion and a Method. The Expansion describes the physical problem, and how the underlying potential is approximated. Currently, DASHMM includes expansions for the Laplace, Yukawa and low-frequency Helmholtz potentials. The Method describes how the approximations to the potential are used to compute the potential at the target locations. DASHMM includes four Methods: the Barnes-Hut method, the Fast Multipole Method, an advanced version of the Fast Multipole Method that uses the merge-and-shift technique (referred to as FMM97), and a direct summation method useful for accuracy comparisons. Not all Expansions are good with every Method. In their current form, each Method in DASHMM has an expansion that is more suited to it.

The next part of `perform_evaluation_test` covers the various pairs of expansion and method:

```C++
  ...

  //Perform the evaluation
  double t0{};
  double tf{};
  int err{0};

  if (args.kernel == std::string{"laplace"}) {
    if (args.method == std::string{"bh"}) {
      dashmm::BH<SourceData, TargetData, dashmm::LaplaceCOM> method{0.6};

      t0 = getticks();
      err = laplace_bh.evaluate(source_handle, target_handle,
                                args.refinement_limit, method,
                                args.accuracy, std::vector<double>{});
      assert(err == dashmm::kSuccess);
      tf = getticks();
    } else if (args.method == std::string{"fmm"}) {
      dashmm::FMM<SourceData, TargetData, dashmm::Laplace> method{};

      t0 = getticks();
      err = laplace_fmm.evaluate(source_handle, target_handle,
                                 args.refinement_limit, method,
                                 args.accuracy, std::vector<double>{});
      assert(err == dashmm::kSuccess);
      tf = getticks();
    } else if (args.method == std::string{"fmm97"}) {
      dashmm::FMM97<SourceData, TargetData, dashmm::Laplace> method{};

      t0 = getticks();
      err = laplace_fmm97.evaluate(source_handle, target_handle,
                                   args.refinement_limit, method,
                                   args.accuracy, std::vector<double>{});
      assert(err == dashmm::kSuccess);
      tf = getticks();
    }
  } else if (args.kernel == std::string{"yukawa"}) {
    if (args.method == std::string{"fmm97"}) {
      dashmm::FMM97<SourceData, TargetData, dashmm::Yukawa> method{};
      std::vector<double> kernelparms(1, 0.1);

      t0 = getticks();
      err = yukawa_fmm97.evaluate(source_handle, target_handle,
                                  args.refinement_limit, method,
                                  args.accuracy, kernelparms);
      assert(err == dashmm::kSuccess);
      tf = getticks();
    }
  } else if (args.kernel == std::string{"helmholtz"}) {
    if (args.method == std::string{"fmm97"}) {
      dashmm::FMM97<SourceData, TargetData, dashmm::Helmholtz> method{};
      std::vector<double> kernelparms(1, 0.1);

      t0 = getticks();
      err = helmholtz_fmm97.evaluate(source_handle, target_handle,
                                     args.refinement_limit, method,
                                     args.accuracy, kernelparms);
      assert(err == dashmm::kSuccess);
      tf = getticks();
    }
  }

  fprintf(stdout, "Evaluation took %lg [us]\n", elapsed(tf, t0));

  ...
```

Each branch above is very similar: the `evaluate()` method is called on one of the `Evaluator` objects defined in this program (see Below).

In DASHMM, Method objects might take parameters during their construction. So `evaluate()` requires an instance of the Method type selected. In the case of FMM and its variation, and the Direct method, these take no arguments. For the Barnes-Hut method, a critical angle must be specified.

For each expansion, there are also two parameters that control the behavior of the expansion provided to `evaluate()`. The first is the accuracy parameter. The meaning of this integer depends on the Expansion being employed, but typically it indicates the number of digits of accuracy requested in the resulting potential computation. The second is a STL vector that contains the kernel parameters for the Expansion in use. These parameters control the exact form of the potential that is being approximated. For details on what is required for each Expansion, please see the documentation.

This just leaves us with the `Evaluator` object, which we turn to next.

## Evaluators

The central object in a DASHMM multipole method evaluation is an `Evaluator` object. These objects provide two major services to the user of DASHMM: they register the needed actions with the underlying HPX-5 runtime system, and they perform the evaluation.

The actions taken during a multipole method evaluation depend on up to four types. These are the same four types that have shown up above: the Source, the Target, the Expansion and the Method. Because DASHMM is intended to be general and extensible, the framework must support not only the built-in Methods and Expansions, but those created by advanced users of the library. This leads to a situation where we have to inform the library before use what exact combinations of these four types are to be used. This allows the library to set up the needed actions.

To inform DASHMM of your intended sets of types, one must create an instance of an `Evaluator` parameterized with the needed types. This must occur before `dashmm::init`. In the demo, these evaluators are created as global objects.

```C++
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::LaplaceCOM, dashmm::BH> laplace_bh{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::LaplaceCOM, dashmm::Direct> laplace_direct{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Laplace, dashmm::FMM> laplace_fmm{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Laplace, dashmm::FMM97> laplace_fmm97{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Yukawa, dashmm::Direct> yukawa_direct{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Yukawa, dashmm::FMM97> yukawa_fmm97{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Helmholtz, dashmm::Direct> helmholtz_direct{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Helmholtz, dashmm::FMM97> helmholtz_fmm97{};
```

The constructor for an `Evaluator` takes no arguments, but the object is a template over four parameters. It is an error to create more than one instance of an `Evaluator` with the same set of four types. There is no state associated with an Evaluator, so this is not a terrible restriction to impose.

The actual evaluation is done through the `evaluate` member of an `Evaluator` object. This call requires six arguments: the `Array` of source data, the `Array` of target data, the refinement limit, an instance of the Method type being employed, an accuracy parameter and a STL vector of kernel parameters. The refinement limit controls the granularity of the hierarchical space partition computed with the method. The source and target points are organized into trees with leaves containing at most a number of particles equal to the refinement limit. For example, a refinement limit of 1 would place a single source or target in the leaves of the tree. Typically, a few tens of points is a good choice, but it is worth experimenting for a specific problem what is a good limit.

## Comparing Accuracy

The final portions of `perform_evaluation_test` run a comparison case with a direct summation method to measure the accuracy of the approximate method.

```C++
  ...

  if (args.verify) {
    // Save a few targets for the direct comparison
    int test_count{0};
    if (hpx_get_my_rank() == 0) {
      test_count = 400;
    }
    if (test_count > args.target_count) {
      test_count = args.target_count;
    }

    //Get the results from the global address space
    args.target_count = target_handle.length();
    TargetData *targets = target_handle.collect();
    if (targets) {
      std::sort(targets, &targets[args.target_count],
                [] (const TargetData &a, const TargetData &b) -> bool {
                  return (a.index < b.index);
                });
    }

    // Copy the test particles into test_targets
    TargetData *test_targets{nullptr};
    if (test_count) {
      test_targets = new TargetData[test_count];

      for (int i = 0; i < test_count; ++i) {
        int idx = i * (args.target_count / test_count);
        assert(idx < args.target_count);
        test_targets[i] = targets[idx];
        test_targets[i].phi = std::complex<double>{0.0, 0.0};
      }
    }

  ...
```

If the user has selected to compare with the exact result computed via the Direct summation method, the first step is to select a small set of points from the targets. This attempts to compare to 400 points, or the full number of target points, whichever is smaller.

During the previous evaluation, the target particles will have been sorted and moved from across ranks in the system. In order to consistently subsample the same target locations, the full set of records is collected onto rank zero. This is performed with the `collect()` method of `Array`. This routine will return a local allocation with all of the records of the distributed array. Once this is accomplished, they are sorted according to their initial index.

After, the sub-sampled target locations are copied into a shorter array which is then copied into a new global `Array` using the now-familiar `allocate()` and `put()` methods:

```C++
  ...

    // Create array for test targets
    dashmm::Array<TargetData> test_handle{};
    err = test_handle.allocate(test_count);
    assert(err == dashmm::kSuccess);
    err = test_handle.put(0, test_count, test_targets);
    assert(err == dashmm::kSuccess);
    delete [] test_targets;

    //do direct evaluation
    if (args.kernel == "laplace") {
      dashmm::Direct<SourceData, TargetData, dashmm::LaplaceCOM> direct{};
      err = laplace_direct.evaluate(source_handle, test_handle,
                                    args.refinement_limit, direct,
                                    args.accuracy, std::vector<double>{});
      assert(err == dashmm::kSuccess);
    } else if (args.kernel == "yukawa") {
      dashmm::Direct<SourceData, TargetData, dashmm::Yukawa> direct{};
      std::vector<double> kernelparms(1, 0.1);
      err = yukawa_direct.evaluate(source_handle, test_handle,
                                   args.refinement_limit, direct,
                                   args.accuracy, kernelparms);
      assert(err == dashmm::kSuccess);
    } else if (args.kernel == "helmholtz") {
      dashmm::Direct<SourceData, TargetData, dashmm::Helmholtz> direct{};
      std::vector<double> kernelparms(1, 0.1);
      err = helmholtz_direct.evaluate(source_handle, test_handle,
                                      args.refinement_limit, direct,
                                      args.accuracy, kernelparms);
      assert(err == dashmm::kSuccess);
    }

    // Retrieve the test results
    test_targets = test_handle.collect();

    //Test error
    compare_results(targets, args.target_count, test_targets, test_count);

    err = test_handle.destroy();
    assert(err == dashmm::kSuccess);
    delete [] test_targets;
    delete [] targets;
  }

  //free up resources
  err = source_handle.destroy();
  assert(err == dashmm::kSuccess);
  err = target_handle.destroy();
  assert(err == dashmm::kSuccess);
}
```

Following the creation of the `test_handle` array, the evaluation is performed, either for the Laplace or Yukawa kernel for the original sources, and the sub-sampled targets. Once again, the results are `collect()`-ed. These are then fed into an error comparison routine `compare_results`, which will not be covered.

## Next Steps

Now you are familiar with a good portion of the elements of the basic interface to DASHMM. Specific details about the return codes from these functions can be found in the DASHMM User Guide, in [Resources]({{ site.baseurl }}/resources.html).

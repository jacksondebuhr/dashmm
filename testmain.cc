#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <map>
#include <memory>

#include <hpx/hpx.h>

#include "dashmm.h"
#include "include/array.h"
#include "include/ids.h"
#include "include/laplace_com.h"
#include "include/reductionops.h"
#include "include/particle.h"


int test_init_handler(int UNUSED) {
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0, test_init_action, test_init_handler, HPX_INT);


int test_array_contents_handler(hpx_addr_t addx) {
  dashmm::ArrayMetaData *meta{nullptr};
  assert(hpx_gas_try_pin(addx, (void **)&meta));
  fprintf(stdout, "Array has %lu entries of size %lu\n",
          meta->count, meta->size);

  double *values{nullptr};
  assert(hpx_gas_try_pin(meta->data, (void **)&values));

  for (size_t i = 0; i < meta->count; ++i) {
    fprintf(stdout, "%lg ", values[i]);
  }
  fprintf(stdout, "\n\n");

  std::sort(values, &values[10]);

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(addx);

  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0,
           test_array_contents_action, test_array_contents_handler,
           HPX_ADDR);


void perform_array_testing(void) {
  dashmm::ObjectHandle array_handle;
  double some_data[10] = {10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};

  auto err = dashmm::allocate_array(10, sizeof(double), &array_handle);
  if (err == dashmm::kRuntimeError) {
    fprintf(stderr, "There was a runtime error during allocate_array\n");
  } else if (err == dashmm::kAllocationError) {
    fprintf(stderr, "There was trouble allocating the requested memory\n");
  } else {
    fprintf(stdout, "Allocation successful\n");
  }

  err = dashmm::array_put(array_handle, 0, 10, some_data);

  //call an action to print that stuff out
  hpx_run(&test_array_contents_action, &array_handle);

  //now pull the data back out
  double sorted_data[10];
  err = dashmm::array_get(array_handle, 0, 10, sorted_data);
  fprintf(stdout, "Pulling the data back out of the global address space:\n");
  for (int i = 0; i < 10; ++i) {
    fprintf(stdout, "%lg ", sorted_data[i]);
  }
  fprintf(stdout, "\n\n");

  err = dashmm::deallocate_array(array_handle);
  if (err == dashmm::kRuntimeError) {
    fprintf(stderr, "There was some trouble in deallocation\n");
  } else {
    fprintf(stdout, "Deallocation was a success!\n");
  }
}


void perform_builtins_testing() {
  //Tests for BH
  dashmm::Method *bhtest = dashmm::bh_method(0.6);
  int bhtype = bhtest->type();
  assert(bhtype == dashmm::kMethodBH);
  //test the release mechanism
  dashmm::MethodSerial *bhserial = bhtest->release();
  assert(bhserial->type == dashmm::kMethodBH);
  assert(bhserial->size == sizeof(double));
  assert(bhserial->data[0] == 0.6);
  delete bhtest;
  //now recreate the method from the data
  auto remade = dashmm::create_method(bhtype, bhserial);
  auto willbenull =
      dashmm::create_method(dashmm::kFirstUserMethodType, bhserial);
  delete bhserial;
  assert(remade->type() == dashmm::kMethodBH);
  assert(willbenull == nullptr);

  {
    dashmm::Expansion *lapcom = dashmm::laplace_COM_expansion();
    dashmm::LaplaceCOM *specific = dynamic_cast<dashmm::LaplaceCOM *>(lapcom);
    specific->set_mtot(1.0);
    specific = nullptr;
    int lapcomtype = lapcom->type();
    size_t lapcombytes = lapcom->bytes();
    assert(lapcomtype == dashmm::kExpansionLaplaceCOM);
    dashmm::LaplaceCOMData *lapcomserial =
        static_cast<dashmm::LaplaceCOMData *>(lapcom->release());
    assert(!lapcom->valid());
    assert(lapcomserial->type == dashmm::kExpansionLaplaceCOM);
    assert(lapcomserial->mtot == 1.0);
    lapcomserial->mtot = 2.0;
    delete lapcom;
    auto remade =
        dashmm::interpret_expansion(lapcomtype, lapcomserial, lapcombytes);
    specific = dynamic_cast<dashmm::LaplaceCOM *>(remade.get());
    assert(specific->term(0).real() == 2.0);
    auto empty =
        dashmm::interpret_expansion(lapcomtype, nullptr, 0);
    auto nowayitworks =
        dashmm::interpret_expansion(dashmm::kFirstUserExpansionType, nullptr, 0);
    assert(remade->type() == lapcomtype);
    assert(remade->bytes() == lapcombytes);
    assert(empty->type() == lapcomtype);
    assert(empty->bytes() == lapcombytes);
    assert(nowayitworks == nullptr);
    delete lapcomserial;

    auto anotherone =
        dashmm::create_expansion(lapcomtype, dashmm::Point{1.0, 2.0, 3.0});
    assert(anotherone->type() == lapcomtype);
    assert(anotherone->bytes() == lapcombytes);
  }
}


int reduction_test_handler(int number) {
  hpx_addr_t red = hpx_lco_reduce_new(number, 2 * sizeof(int),
                                dashmm::int_sum_ident_op, dashmm::int_sum_op);
  for (int i = 1; i <= number; ++i) {
    int vals[2] = {i, 2*i};
    hpx_lco_set_lsync(red, 2 * sizeof(int), vals, HPX_NULL);
  }
  int results[2];
  hpx_lco_get(red, 2 * sizeof(int), results);
  int compare = (number * (number + 1)) / 2;
  assert(results[0] == compare);
  assert(results[1] == compare * 2);

  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0,
           reduction_test_action, reduction_test_handler, HPX_INT);


void perform_reduction_testing() {
  int terms = 10;
  hpx_run(&reduction_test_action, &terms);
}


int particle_testing_handler(int UNUSED) {
  //both the source and target refs...

  //create some sources
  dashmm::Source src[2] = {{1.0, dashmm::Point{0.0, 0.0, 0.0}},
                           {1.0, dashmm::Point{1.0, 0.0, 0.0}}};
  dashmm::SourceRef srcref{src, 2};
  dashmm::SourceRef otherref{srcref.data(), srcref.n()};

  //create some targets
  dashmm::Target trg[2] = {
      {dashmm::Point{0.5, 0.0, 0.0}, std::complex<double>(0.0, 0.0), 0},
      {dashmm::Point{0.0, 1.0, 0.0}, std::complex<double>(0.0, 0.0), 1}};
  dashmm::TargetRef trgref{trg, 2};
  dashmm::TargetRef anotherref{trgref.data(), trgref.n()};
  trgref.schedule(2);
  trgref.finalize();
  trgref.contribute_S_to_T(dashmm::kExpansionLaplaceCOM, 2, src);
  dashmm::LaplaceCOM lapcom{dashmm::Point{0.0, 0.0, 0.0}};
  lapcom.set_mtot(2.0);
  double xcom[3] = {0.5, 1.0, 0.0};
  double Q[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  lapcom.set_xcom(xcom);
  lapcom.set_Q(Q);
  dashmm::LaplaceCOMData *lapcomserial =
      (dashmm::LaplaceCOMData *)lapcom.release();
  anotherref.contribute_M_to_T(lapcom.type(), lapcom.bytes(), lapcomserial);
  free(lapcomserial);
  hpx_lco_wait(anotherref.data());
  struct TargetRefLCOData {
    int arrived;
    int scheduled;
    int finished;
    int count;
    dashmm::Target targets[];
  } *tdata;
  hpx_lco_getref(anotherref.data(), 1, (void **)&tdata);
  //check results
  assert(tdata->arrived == 2);
  assert(tdata->scheduled == 2);
  assert(tdata->finished == 1);
  assert(tdata->count == 2);
  fprintf(stdout, "Target 0: (%lg %lg %lg) - %lg\n",
          tdata->targets[0].position.x(), tdata->targets[0].position.y(),
          tdata->targets[0].position.z(), tdata->targets[0].phi.real());
  fprintf(stdout, "Target 1: (%lg %lg %lg) - %lg\n",
          tdata->targets[1].position.x(), tdata->targets[1].position.y(),
          tdata->targets[1].position.z(), tdata->targets[1].phi.real());
  hpx_lco_release(trgref.data(), tdata);

  //clean up
  srcref.destroy();
  trgref.destroy();

  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0,
           particle_testing_action, particle_testing_handler,
           HPX_INT);


void perform_particle_testing() {
  int unused = 42;
  hpx_run(&particle_testing_action, &unused);
}


constexpr int source_test_count = 10000;
constexpr int target_test_count = 10000;

struct UserSourceData {
  double pos[3];
  double mass;
};

struct UserTargetData {
  double pos[3];
  double phi[2];    //real, imag
  double phi_direct[2];
};


void perform_the_big_test() {
  //create some arrays
  UserSourceData *sources = static_cast<UserSourceData *>(
        malloc(sizeof(UserSourceData) * source_test_count));
  UserTargetData *targets = static_cast<UserTargetData *>(
        malloc(sizeof(UserTargetData) * target_test_count));
  srand(12345);
  for (int i = 0; i < source_test_count; ++i) {
    sources[i].pos[0] = (double)rand() / RAND_MAX;
    sources[i].pos[1] = (double)rand() / RAND_MAX;
    sources[i].pos[2] = (double)rand() / RAND_MAX;
    sources[i].mass = (double)rand() / RAND_MAX + 1.0;
  }
  for (int i = 0; i < target_test_count; ++i) {
    targets[i].pos[0] = (double)rand() / RAND_MAX;
    targets[i].pos[1] = (double)rand() / RAND_MAX;
    targets[i].pos[2] = (double)rand() / RAND_MAX;
    targets[i].phi[0] = 0.0;
    targets[i].phi[1] = 0.0;
    targets[i].phi_direct[0] = 0.0;
    targets[i].phi_direct[1] = 0.0;
  }

  //prep sources
  dashmm::ObjectHandle source_handle;
  auto err = dashmm::allocate_array(source_test_count, sizeof(UserSourceData),
                            &source_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::array_put(source_handle, 0, source_test_count, sources);
  assert(err == dashmm::kSuccess);

  //prep targets
  dashmm::ObjectHandle target_handle;
  err = dashmm::allocate_array(target_test_count, sizeof(UserTargetData),
                               &target_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::array_put(target_handle, 0, target_test_count, targets);
  assert(err == dashmm::kSuccess);

  //get method and expansion
  auto bhmethod = dashmm::bh_method(0.6);
  assert(bhmethod != nullptr);
  auto lapcomexp = dashmm::laplace_COM_expansion();
  assert(lapcomexp != nullptr);

  //evaluate
  //TODO: Is it weird that evaluate eats the method and expansion?
  //*
  err = dashmm::evaluate(source_handle, offsetof(UserSourceData, pos),
                         offsetof(UserSourceData, mass),
                         target_handle, offsetof(UserTargetData, pos),
                         offsetof(UserTargetData, phi),
                         40,
                         std::unique_ptr<dashmm::Method>{bhmethod},
                         std::unique_ptr<dashmm::Expansion>{lapcomexp});
  //*/

  //*
  auto direct = dashmm::bh_method(0.0001);
  auto direxp = dashmm::laplace_COM_expansion();
  err = dashmm::evaluate(source_handle, offsetof(UserSourceData, pos),
                         offsetof(UserSourceData, mass),
                         target_handle, offsetof(UserTargetData, pos),
                         offsetof(UserTargetData, phi_direct),
                         40,
                         std::unique_ptr<dashmm::Method>{direct},
                         std::unique_ptr<dashmm::Expansion>{direxp});
  //*/

  //get targets
  err = dashmm::array_get(target_handle, 0, target_test_count, targets);
  assert(err == dashmm::kSuccess);
  //TODO: check that the data comes out
  for (int i = 0; i < target_test_count; ++i) {
    fprintf(stdout, "%d: %lg %lg %lg\n", i,
        targets[i].phi[0], targets[i].phi_direct[0],
        fabs(targets[i].phi[0] - targets[i].phi_direct[0])
                / targets[i].phi_direct[0]);
  }

  //free up stuff
  err = dashmm::deallocate_array(source_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::deallocate_array(target_handle);
  assert(err == dashmm::kSuccess);

  free(sources);
  free(targets);
}


int main(int argc, char **argv) {
  int perform_other_tests = 0;

  auto err = dashmm::init(&argc, &argv);
  if (err == dashmm::kInitError) {
    fprintf(stderr, "dashmm::init gave an error!\n");
  } else if (err == dashmm::kRuntimeError) {
    fprintf(stderr, "dashmm::init had trouble starting runtime\n");
  } else if (err != dashmm::kSuccess) {
    fprintf(stderr, "Woah! init() gave an unknown error\n");
  }

  //Do a quick test that init did something.
  if (perform_other_tests) {
    int throwaway = 42;
    hpx_run(&test_init_action, &throwaway);

    //now let's work on the array stuff
    perform_array_testing();

    //builtins testing
    perform_builtins_testing();

    //test reduction ops
    perform_reduction_testing();

    //test the particle types
    perform_particle_testing();
  }


  perform_the_big_test();


  err = dashmm::finalize();
  if (err == dashmm::kFiniError) {
    fprintf(stderr, "dashmm::finalize gave an error\n");
  } else if (err != dashmm::kSuccess) {
    fprintf(stderr, "Woah! finalize() gave an unknown error\n");
  }

  return 0;
}

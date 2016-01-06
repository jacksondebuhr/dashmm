#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <map>
#include <memory>

#include "dashmm.h"


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


  perform_the_big_test();


  err = dashmm::finalize();
  if (err == dashmm::kFiniError) {
    fprintf(stderr, "dashmm::finalize gave an error\n");
  } else if (err != dashmm::kSuccess) {
    fprintf(stderr, "Woah! finalize() gave an unknown error\n");
  }

  return 0;
}

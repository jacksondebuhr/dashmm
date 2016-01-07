#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <algorithm>
#include <map>
#include <memory>

#include "dashmm.h"


struct UserSourceData {
  double pos[3];
  double mass;
};

struct UserTargetData {
  double pos[3];
  double phi[2];    //real, imag
};


void pick_cube_position(double *pos) {
  pos[0] = (double)rand() / RAND_MAX;
  pos[1] = (double)rand() / RAND_MAX;
  pos[2] = (double)rand() / RAND_MAX;
}


void pick_sphere_position(double *pos) {
  double r = 1.0;
  double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0;
  double stheta = sqrt(1.0 - ctheta * ctheta);
  double phi = 2.0 * 3.1415926535 * (double)rand() / RAND_MAX;
  pos[0] = r * stheta * cos(phi);
  pos[1] = r * stheta * sin(phi);
  pos[2] = r * ctheta;
}


double pick_mass() {
  return (double)rand() / RAND_MAX + 1.0;
}


void pick_plummer_position(double *pos) {
  //NOTE: This is using a = 1
  double unif = (double)rand() / RAND_MAX;
  double r = 1.0 / sqrt(pow(unif, -2.0 / 3.0) - 1);
  double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0;
  double stheta = sqrt(1.0 - ctheta * ctheta);
  double phi = 2.0 * 3.1415926535 * (double)rand() / RAND_MAX;
  pos[0] = r * stheta * cos(phi);
  pos[1] = r * stheta * sin(phi);
  pos[2] = r * ctheta;
}


double pick_plummer_mass(int count) {
  //We take the total mass to be 100
  return 100.0 / count;
}


void set_sources(UserSourceData *sources, int source_count, int source_type) {
  if (source_type == 0) {
    //Cube
    for (int i = 0; i < source_count; ++i) {
      pick_cube_position(sources[i].pos);
      sources[i].mass = pick_mass();
    }
  } else if (source_type == 1) {
    //Sphere
    for (int i = 0; i < source_count; ++i) {
      pick_sphere_position(sources[i].pos);
      sources[i].mass = pick_mass();
    }
  } else {
    //Plummer
    for (int i = 0; i < source_count; ++i) {
      pick_plummer_position(sources[i].pos);
      sources[i].mass = pick_plummer_mass(source_count);
    }
  }
}


void set_targets(UserTargetData *targets, int target_count, int target_type) {
  if (target_type == 0) {
    //Cube
    for (int i = 0; i < target_count; ++i) {
      pick_cube_position(targets[i].pos);
      targets[i].phi[0] = 0.0;
      targets[i].phi[1] = 0.0;
    }
  } else if (target_type == 1) {
    //Sphere
    for (int i = 0; i < target_count; ++i) {
      pick_sphere_position(targets[i].pos);
      targets[i].phi[0] = 0.0;
      targets[i].phi[1] = 0.0;
    }
  } else {
    //Plummer
    for (int i = 0; i < target_count; ++i) {
      pick_plummer_position(targets[i].pos);
      targets[i].phi[0] = 0.0;
      targets[i].phi[1] = 0.0;
    }
  }
}


void perform_evaluation_test(int source_count, int source_type,
                             int target_count, int target_type,
                             int refinement_limit, int test_case) {
  srand(123456);

  //work out the test count
  int test_count = 400;
  if (test_count > target_count) {
    test_count = target_count;
  }

  //create some arrays
  UserSourceData *sources = static_cast<UserSourceData *>(
        malloc(sizeof(UserSourceData) * source_count));
  UserTargetData *targets = static_cast<UserTargetData *>(
        malloc(sizeof(UserTargetData) * target_count));
  UserTargetData *test_targets = static_cast<UserTargetData *>(
        malloc(sizeof(UserTargetData) * test_count));

  set_sources(sources, source_count, source_type);
  set_targets(targets, target_count, target_type);
  //Save a few targets for direct comparison
  for (int i = 0; i < test_count; ++i) {
    int idx = i * (target_count / test_count);
    assert(idx < target_count);
    test_targets[i] = targets[idx];
  }

  //prep sources
  dashmm::ObjectHandle source_handle;
  auto err = dashmm::allocate_array(source_count, sizeof(UserSourceData),
                            &source_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::array_put(source_handle, 0, source_count, sources);
  assert(err == dashmm::kSuccess);

  //prep targets
  dashmm::ObjectHandle target_handle;
  err = dashmm::allocate_array(target_count, sizeof(UserTargetData),
                               &target_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::array_put(target_handle, 0, target_count, targets);
  assert(err == dashmm::kSuccess);

  //prep test targets
  dashmm::ObjectHandle test_handle;
  err = dashmm::allocate_array(test_count, sizeof(UserTargetData),
                               &test_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::array_put(test_handle, 0, test_count, test_targets);
  assert(err == dashmm::kSuccess);

  //get method and expansion
  //TODO: generalize based on the inputs
  dashmm::Method *test_method{nullptr};
  dashmm::Expansion *test_expansion{nullptr};
  if (test_case == 0) {
    test_method = dashmm::bh_method(0.6);
    test_expansion = dashmm::laplace_COM_expansion();
  } //more cases once they are implemented
  assert(test_method && test_expansion);

  //evaluate - first for the approximate version
  err = dashmm::evaluate(source_handle, offsetof(UserSourceData, pos),
                         offsetof(UserSourceData, mass),
                         target_handle, offsetof(UserTargetData, pos),
                         offsetof(UserTargetData, phi),
                         refinement_limit,
                         std::unique_ptr<dashmm::Method>{test_method},
                         std::unique_ptr<dashmm::Expansion>{test_expansion});

  //This is effectively the exact potential computation - these are not
  // subject to the input case_type because we just want to do the direct
  // computation here.
  //TODO: change this to use the direct summation method, once that is
  // implemented
  auto direct = dashmm::bh_method(0.00001);
  auto direxp = dashmm::laplace_COM_expansion();
  err = dashmm::evaluate(source_handle, offsetof(UserSourceData, pos),
                         offsetof(UserSourceData, mass),
                         test_handle, offsetof(UserTargetData, pos),
                         offsetof(UserTargetData, phi),
                         refinement_limit,
                         std::unique_ptr<dashmm::Method>{direct},
                         std::unique_ptr<dashmm::Expansion>{direxp});

  //get targets
  err = dashmm::array_get(target_handle, 0, target_count, targets);
  assert(err == dashmm::kSuccess);
  err = dashmm::array_get(test_handle, 0, test_count, test_targets);
  assert(err == dashmm::kSuccess);

  //Test error
  double numerator = 0.0;
  double denominator = 0.0;
  double maxrel = 0.0;
  for (int i = 0; i < test_count; ++i) {
    int idx = i * (target_count / test_count);
    assert(idx < target_count);
    double relerr = fabs(targets[idx].phi[0] - test_targets[i].phi[0]);
    numerator += relerr * relerr;
    denominator += test_targets[i].phi[0] * test_targets[i].phi[0];
    if (relerr / test_targets[i].phi[0] > maxrel) {
      maxrel = relerr / test_targets[i].phi[0];
    }
  }
  fprintf(stdout, "Error for %d test points: %lg (max %lg)\n",
                  test_count, sqrt(numerator / denominator), maxrel);

  //free up stuff
  err = dashmm::deallocate_array(source_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::deallocate_array(target_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::deallocate_array(test_handle);
  assert(err == dashmm::kSuccess);

  free(sources);
  free(targets);
  free(test_targets);
}


int main(int argc, char **argv) {
  auto err = dashmm::init(&argc, &argv);
  if (err == dashmm::kInitError) {
    fprintf(stderr, "dashmm::init gave an error!\n");
    return -1;
  } else if (err == dashmm::kRuntimeError) {
    fprintf(stderr, "dashmm::init had trouble starting runtime\n");
    return -1;
  } else if (err != dashmm::kSuccess) {
    fprintf(stderr, "Woah! init() gave an unknown error\n");
    return -1;
  }


  //TODO: test usage


  //TODO: collect arguments


  perform_evaluation_test(10000, 0, 10000, 0, 40, 0);


  err = dashmm::finalize();
  if (err == dashmm::kFiniError) {
    fprintf(stderr, "dashmm::finalize gave an error\n");
  } else if (err != dashmm::kSuccess) {
    fprintf(stderr, "Woah! finalize() gave an unknown error\n");
  }

  return 0;
}

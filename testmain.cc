#include <cstdio>

#include <algorithm>
#include <map>

#include "dashmm.h"
#include "include/testing.h"
#include "include/array.h"


int test_init_handler(int UNUSED) {
  dashmm::print_method_table();
  dashmm::print_expansion_table();
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
  double some_data[100] = {10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};

  auto err = dashmm::allocate_array(10, sizeof(double), &array_handle);
  if (err == dashmm::kRuntimeError) {
    fprintf(stderr, "There was a runtime error during allocate_array\n");
  } else if (err == dashmm::kAllocationError) {
    fprintf(stderr, "There was trouble allocating the requested memory\n");
  } else {
    fprintf(stdout, "Allocation successful\n");
  }

  err = dashmm::array_put(array_handle, 0, 20, some_data);

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


int main(int argc, char **argv) {
  auto err = dashmm::init(&argc, &argv);
  if (err == dashmm::kInitError) {
    fprintf(stderr, "dashmm::init gave an error!\n");
  } else if (err == dashmm::kRuntimeError) {
    fprintf(stderr, "dashmm::init had trouble starting runtime\n");
  } else if (err != dashmm::kSuccess) {
    fprintf(stderr, "Woah! init() gave an unknown error\n");
  }


  //Do a quick test that init did something.
  int throwaway = 42;
  hpx_run(&test_init_action, &throwaway);


  //now let's work on the array stuff
  perform_array_testing();


  err = dashmm::finalize();
  if (err == dashmm::kFiniError) {
    fprintf(stderr, "dashmm::finalize gave an error\n");
  } else if (err != dashmm::kSuccess) {
    fprintf(stderr, "Woah! finalize() gave an unknown error\n");
  }

  return 0;
}

#include <iostream>

#include <hpx/hpx.h>


template <typename T>
int accumulate(T *vec, size_t bytes) {
  size_t count = bytes / sizeof(T);
  std::cout << "I receive: ";
  T total{};
  for (size_t i = 0; i < count; ++i) {
    total += vec[i];
    std::cout << vec[i] << " ";
  }
  std::cout << std::endl << "I compute: " << total << std::endl;
  hpx_exit(HPX_SUCCESS);
}


HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           accumulate_double_action, accumulate<double>,
           HPX_POINTER, HPX_SIZE_T);

HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           accumulate_int_action, accumulate<int>,
           HPX_POINTER, HPX_SIZE_T);


int main(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) {
    return -1;
  }

  //call template action
  double args[3] = {1.0, 2.2, 3.7};
  hpx_run(&accumulate_double_action, args, sizeof(double) * 3);

  int argsi[6] = {1, 2, 3, 4, 5, 6};
  hpx_run(&accumulate_int_action, argsi, sizeof(int) * 6);

  hpx_finalize();
  return 0;
}

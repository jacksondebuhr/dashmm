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


#define SOME_THING(type)                                        \
  HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,                       \
           accumulate_##type##_action, accumulate<type>,        \
           HPX_POINTER, HPX_SIZE_T)                             \

#define ACTION_NAME(type) accumulate_##type##_action

template <typename T>
struct better_action_name {
  const hpx_addr_t value = ACTION_NAME(T);
};


SOME_THING(int);
SOME_THING(double);


template <typename T>
void did_it_work(T *args, size_t bytes) {
  hpx_run(&better_action_name<T>::value, args, bytes);
}


int main(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) {
    return -1;
  }

  //call template action
  double args[3] = {1.0, 2.2, 3.7};
  //hpx_run(&ACTION_NAME(double), args, sizeof(double) * 3);
  did_it_work(args, sizeof(double) * 3);

  int argsi[6] = {1, 2, 3, 4, 5, 6};
  hpx_run(&ACTION_NAME(int), argsi, sizeof(int) * 6);

  hpx_finalize();
  return 0;
}

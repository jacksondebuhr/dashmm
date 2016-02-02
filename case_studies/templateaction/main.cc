#include <iostream>

#include <hpx/hpx.h>


template <typename T, template <typename> class Op>
class dashmm {
 private:
   hpx_action_t the_action_id;

   static int accumulate_handler(T *vec, size_t bytes) {
     size_t count = bytes / sizeof(T);
     std::cout << "I receive: ";
     T total = Op<T>::initial_value();
     for (size_t i = 0; i < count; ++i) {
       Op<T>::reduce(total, vec[i]);
       std::cout << vec[i] << " ";
     }
     std::cout << std::endl << "I compute: " << total << std::endl;
     hpx_exit(HPX_SUCCESS);
   }

 public:
   dashmm() {
     HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, the_action_id,
                         accumulate_handler, HPX_POINTER, HPX_SIZE_T);
   }

  int accumulate(T *vec, size_t num) {
    return hpx_run(&the_action_id, vec, sizeof(T) * num);
  }
};


template <typename T>
struct SumOp {
  static void reduce(T &total, T term) {
    total += term;
  }
  static T initial_value() {
    return 0;
  }
};

template <typename T>
struct MultOp {
  static void reduce(T &total, T term) {
    total *= term;
  }
  static T initial_value() {
    return 1;
  }
};


dashmm<int, SumOp> dashforint{};
dashmm<double, MultOp> dashfordouble{};


int main(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) {
    return -1;
  }

  int argsi[6] = {1, 2, 3, 4, 5, 6};
  dashforint.accumulate(argsi, 6);

  double args[4] = {2.0, 3.14, 1.2245, 11.0};
  dashfordouble.accumulate(args, 4);

  hpx_finalize();
  return 0;
}

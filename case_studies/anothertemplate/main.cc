#include <iostream>

#include <hpx/hpx.h>

#include "evaluator.h"


Evaluator<int> dashforint{};
Evaluator<double> dashfordouble{};


int main(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) {
    return -1;
  }

  int argsi[2] = {10, 15};
  dashforint.accumulate(argsi, 2);

  double args[2] = {2.0, 3.14};
  dashfordouble.accumulate(args, 2);

  hpx_finalize();
  return 0;
}

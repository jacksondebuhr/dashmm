#include "include/reductionops.h"


namespace dashmm {


void int_sum_ident_handler(int *input, const size_t bytes) {
  int count = bytes / sizeof(int);
  for (int i = 0; i < count; ++i) {
    input[i] = 0;
  }
}
HPX_ACTION(HPX_FUNCTION, 0, int_sum_ident_op, int_sum_ident_handler,
           HPX_POINTER, HPX_SIZE_T);


void int_sum_op_handler(int *lhs, const int *rhs, size_t bytes) {
  int count = bytes / sizeof(int);
  for (int i = 0; i < count; ++i) {
    lhs[i] += rhs[i];
  }
}
HPX_ACTION(HPX_FUNCTION, 0, int_sum_op, int_sum_op_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


} // namespace dashmm

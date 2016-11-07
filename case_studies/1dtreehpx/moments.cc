#include "moments.h"


void moment_reduction_ident_handler(Moment *ident, size_t bytes) {
  ident->mtot = 0.0;
  ident->xcom = 0.0;
  ident->Q00  = 0.0;
}
HPX_ACTION(HPX_FUNCTION, 0,
           moment_reduction_ident, moment_reduction_ident_handler,
           HPX_POINTER, HPX_SIZE_T);


void moment_reduction_op_handler(Moment *lhs,
                                 const Moment *rhs, size_t bytes) {
  double newtot = lhs->mtot + rhs->mtot;
  double newcom{0.0};
  double newQ00{0.0};
  if (newtot > 0.0) {
    newcom = lhs->mtot * lhs->xcom + rhs->mtot * rhs->xcom;
    newcom /= newtot;

    double dxl{lhs->xcom - newcom};
    newQ00 = 2.0 * lhs->mtot * dxl * dxl + lhs->Q00;

    double dxr{rhs->xcom - newcom};
    newQ00 += 2.0 * rhs->mtot * dxr * dxr + rhs->Q00;
  }

  lhs->mtot = newtot;
  lhs->xcom = newcom;
  lhs->Q00  = newQ00;
}
HPX_ACTION(HPX_FUNCTION, 0,
           moment_reduction_op, moment_reduction_op_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
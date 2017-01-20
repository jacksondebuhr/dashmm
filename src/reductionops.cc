// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


/// \file
/// \brief Implemention of common reduction operations


#include <limits>
#include <algorithm>
#include "dashmm/reductionops.h"


namespace dashmm {


void int_sum_ident_handler(int *input, const size_t bytes) {
  size_t count = bytes / sizeof(int);
  for (size_t i = 0; i < count; ++i) {
    input[i] = 0;
  }
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, int_sum_ident_op, int_sum_ident_handler,
           HPX_POINTER, HPX_SIZE_T);


void int_sum_op_handler(int *lhs, const int *rhs, size_t bytes) {
  size_t count = bytes / sizeof(int);
  for (size_t i = 0; i < count; ++i) {
    lhs[i] += rhs[i];
  }
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, int_sum_op, int_sum_op_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


void size_sum_ident_handler(size_t *input, const size_t bytes) {
  int count = bytes / sizeof(size_t);
  for (int i = 0; i < count; ++i) {
    input[i] = 0;
  }
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, size_sum_ident, size_sum_ident_handler,
           HPX_POINTER, HPX_SIZE_T);


void size_sum_op_handler(size_t *lhs, const size_t *rhs, size_t bytes) {
  int count = bytes / sizeof(size_t);
  for (int i = 0; i < count; ++i) {
    lhs[i] += rhs[i];
  }
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, size_sum_op, size_sum_op_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

void int_max_ident_handler(int *input, const size_t bytes) {
  size_t count = bytes / sizeof(int);
  for (size_t i = 0; i < count; ++i) {
    input[i] = std::numeric_limits<int>::min();
  }
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, int_max_ident_op, int_max_ident_handler,
           HPX_POINTER, HPX_SIZE_T);

void int_max_op_handler(int *lhs, const int *rhs, size_t bytes) {
  size_t count = bytes / sizeof(int);
  for (size_t i = 0; i < count; ++i) {
    lhs[i] = std::max(lhs[i], rhs[i]);
  }
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, int_max_op, int_max_op_handler,
           HPX_POINTER, HPX_SIZE_T);


} // namespace dashmm

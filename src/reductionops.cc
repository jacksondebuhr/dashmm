// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


/// \file
/// \brief Implemention of common reduction operations


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


} // namespace dashmm

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


#ifndef __DASHMM_HILBERT_H__
#define __DASHMM_HILBERT_H__


namespace dashmm {


/// Distribute the nodes according to a Hilbert space filling curve
int *distribute_points_hilbert(int num_ranks,
                               const int *global,
                               int len,
                               int lvl);


} // dashmm


#endif
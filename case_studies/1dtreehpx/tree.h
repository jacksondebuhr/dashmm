// ============================================================================
//  High Performance ParalleX Library (hpx-apps)
//
//  Copyright (c) 2013-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license.  See the COPYING file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// ============================================================================

#ifndef __1DTREE_TUTORIAL_TREE_H__
#define __1DTREE_TUTORIAL_TREE_H__


#include "hpx/hpx.h"


void broadcast_domain_size(double domain_size);


hpx_addr_t create_node(double low, double high);

void partition_node_sync(hpx_addr_t node, hpx_addr_t parts, int n_parts,
                         int n_partition);

void spawn_potential_computation(hpx_addr_t root, hpx_addr_t sync,
                                 double theta);


#endif // __1DTREE_TUTORIAL_TREE_H__

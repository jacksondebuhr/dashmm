// =============================================================================
//  DASHMM 
//
//  Draft of tree construction for DASHMM
//  dashmm-tree.h
//
//  Copyright (c) 2014, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license.  See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
//
//  Authors:
//    Jackson DeBuhr, Indiana University <jdebuhr [at] indiana.edu>
// =============================================================================

#ifndef __DASHMM_TREE_H__
#define __DASHMM_TREE_H__


#include "hpx/hpx.h"


///
/// \type dashmm_volume_t
///
/// Used to represent rectangular volumes.
///
typedef struct {
  double a[3];          ///< The low corner of the rectangular prism
  double b[3];          ///< The high corner of the rectangular prism
} dashmm_volume_t;


///
/// \type dashmm_point_t
///
/// The point type will be used to represent the source or target points 
/// for DASHMM. In this incarnation, only the positions and an integer used
/// during sorting are saved.
///
/// Eventually, the moments or the forces would be added to this type.
///
typedef struct {
  double pos[3];        ///< The position of the point
  int sort;             ///< Used in sorting the points geometrically
} dashmm_point_t;


///
/// \type dashmm_tree_node_t
///
/// The nodes of the tree. Most of the members are self-explanatory. In a few
/// cases some members are used for multiple purposes. These cases are
/// documented below. In the future, likely the multiple use cases will be 
/// captured with a union.
///
/// Eventually the nodes will also need to store moment information and likely
/// other values.
/// 
typedef struct {
  dashmm_volume_t vol;  ///< the volume of the node
  
  hpx_addr_t parent;    ///< parent of this node
  hpx_addr_t child[8];  ///< children of this node
  
  hpx_addr_t points;    ///< address of the first point represented by this node
  uint32_t n_points;    ///< number of points represented by this node
} dashmm_tree_node_t;


///
/// \type dashmm_tree_t
///
/// This type stores all the relevant information about a tree. This includes
/// the address of the topnodes, and the various parameters that define
/// the resulting tree.
typedef struct {
  dashmm_volume_t vol;  ///< the volume of the root
  int refinement;       ///< the node refinement limit
  int top_depth;        ///< the depth of the top of the tree
  hpx_addr_t topnodes;  ///< the topnodes of the tree
  int n_topnodes;       ///< the number of topnodes
} dashmm_tree_t;



//
// User Interface
//


///
/// \brief Create a tree from the given points.
///
/// The returned address will point to some global memory containing the
/// data in a dashmm_tree_t. This function will not return before the tree is
/// complete, and the data is globally visible in the returned address.
///
/// The tree itself will be in two parts, cut based on depth from the root 
/// node. The top part of the tree, down to @param top_depth, will be a 
/// complete refinement of the given volume. These `topnodes' will always
/// exist in the tree, and could be empty of points. Below those are the
/// `branches' of the tree. Associated with each of the most refined topnodes
/// is a branch of nodes that refine the contined particles to the given
/// @param refinement limit. The incoming point data is replicated and 
/// reordered in the branches of the tree. This memory replication is intended
/// to model the eventual use case that the input points data is taken from
/// a more complete user-defined record that may contain application specific
/// information about which DASHMM does not care.
///
/// However, to work with the input points, the system needs to know the 
/// distribution of those points in the global address space, hence the
/// @param n_per_block parameter.
///
/// \param points - the global address of the points from which to form the tree
/// \param n_points - the number of points
/// \param n_per_block - the number of points per GAS block in @param points
/// \param refinement - the largest number of points allowed in a leaf node
/// \param top_depth - the number of levels of full refinement at the top of
///                    the tree.
/// 
/// \return - The global address of the created tree; the tree will be in the
///           global address space, but the calling thread will need to share
///           the address if other threads need the address
///
hpx_addr_t dashmm_tree_create(hpx_addr_t points,
                              int n_points,
                              int n_per_block,
                              int refinement,
                              int top_depth);
                              

///
/// \brief Destroy a dashmm tree
///
/// This will destroy the tree, and all memory associated with it. This routine
/// proceeds immediately, and so it is imperative that the tree not be in use
/// when this routine is called.
///
/// This will remove all branches of the tree and the point data allocated to
/// store the sorted list of points.
///
/// \param tree_gas - the global address of the tree to destroy
///
void dashmm_tree_destroy(hpx_addr_t tree_gas);



//
// Utility functions
//


void dashmm_tree_register_actions(void);


int dashmm_volume_which_child(dashmm_volume_t vol, double *pos);
dashmm_volume_t dashmm_volume_of_child(dashmm_volume_t vol, int which);
void dashmm_volume_reducer(void *a, void *b);
void dashmm_volume_cubify(dashmm_volume_t *vol, double expand);


void dashmm_tree_topnode_index_in_level(int index, int level, int *array);
int dashmm_tree_topnode_offset_on_level(int level, int *array);
int dashmm_tree_topnode_offset_from_index(int level, int *array);
dashmm_volume_t dashmm_tree_topnode_volume(dashmm_volume_t vol, 
                                           int level, int *index);
hpx_addr_t dashmm_tree_topnode_address(hpx_addr_t base, int level, int *index);
void dashmm_tree_topnode_from_point(dashmm_volume_t vol, int level,
                                    dashmm_point_t *point, int *index);


void dashmm_point_bin_sort(dashmm_point_t *points, int n_points, 
                            int *bins, int n_bins);


void dashmm_continue_cleanup(void *env);


//NOTES: the caller needs to free the returned future
// This works for actions whose parameters are the block size and block number.
//
// We need to upgrade this to take other parameters as well...
// They should be constant over the blocks. This is needed for the sorting case,
// because we need to have access to the topnode array so that we can target
// the actions correctly.
hpx_addr_t dashmm_parallel_block_spawn(hpx_addr_t base,
                                       int block_size,
                                       int block_count,
                                       hpx_action_t terminal_action,
                                       int payload_size,
                                       void *payload,
                                       int result_size,
                                       void (*reduce)(void *a, void *b) );


//
// Actions
//


int dashmm_block_spawn(void *args);

int dashmm_point_volume(void *args);

int dashmm_tree_node_init(void *args);

int dashmm_tree_points_sort(void *args);

int dashmm_tree_inform_topnode(void *args);

int dashmm_tree_node_alloc(void *args);

int dashmm_tree_node_refine(void *args);

int dashmm_tree_points_refine(void *args);

int dashmm_tree_init_child(void *args);


#endif

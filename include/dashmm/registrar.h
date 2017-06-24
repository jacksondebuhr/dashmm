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


#ifndef __DASHMM_REGISTRAR_H__
#define __DASHMM_REGISTRAR_H__


/// \file
/// \brief Registrar objects for HPX-5 active objects


#include "dashmm/expansionlco.h"
#include "dashmm/targetlco.h"
#include "dashmm/tree.h"


namespace dashmm {


/// Object that handles action registration for Array objects
template <typename T>
class ArrayRegistrar {
 public:
  using array_t = Array<T>;
  ArrayRegistrar() {
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::array_local_count_,
                        array_t::array_local_count_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::array_total_count_,
                        array_t::array_total_count_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::allocate_array_meta_,
                        array_t::allocate_array_meta_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::allocate_local_work_,
                        array_t::allocate_local_work_handler,
                        HPX_ADDR, HPX_ADDR, HPX_ADDR,
                        HPX_SIZE_T, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::allocate_array_destroy_reducer_,
                        array_t::allocate_array_destroy_reducer_handler,
                        HPX_ADDR, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::deallocate_array_,
                        array_t::deallocate_array_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::deallocate_array_local_,
                        array_t::deallocate_array_local_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::get_or_put_retcode_reducer_,
                        array_t::get_or_put_retcode_reducer_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::get_or_put_reducer_delete_,
                        array_t::get_or_put_reducer_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::array_get_,
                        array_t::array_get_handler,
                        HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T,
                        HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::array_put_,
                        array_t::array_put_handler,
                        HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T,
                        HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::array_segment_request_,
                        array_t::array_segment_request_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::array_collect_,
                        array_t::array_collect_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::array_collect_request_,
                        array_t::array_collect_request_handler,
                        HPX_ADDR, HPX_ADDR, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        array_t::array_collect_receive_,
                        array_t::array_collect_receive_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        array_t::array_set_manager_,
                        aaray_t::array_set_manager_handler,
                        HPX_POINTER, HPX_ADDR);
  }
};


/// Object that handles action registration for TargetLCOs
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class TargetLCORegistrar {
 public:
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;

  TargetLCORegistrar() {
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        targetlco_t::init_, targetlco_t::init_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        targetlco_t::operation_,
                        targetlco_t::operation_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        targetlco_t::predicate_,
                        targetlco_t::predicate_handler,
                        HPX_POINTER, HPX_SIZE_T);
  }
};


/// Object that handles action registration for ExpansionLCOs
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class ExpansionLCORegistrar {
 public:
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;

  ExpansionLCORegistrar() {
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        expansionlco_t::init_, expansionlco_t::init_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        expansionlco_t::operation_,
                        expansionlco_t::operation_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        expansionlco_t::predicate_,
                        expansionlco_t::predicate_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        expansionlco_t::spawn_out_edges_,
                        expansionlco_t::spawn_out_edges_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::spawn_out_edges_from_remote_,
                        expansionlco_t::spawn_out_edges_from_remote_handler,
                        HPX_POINTER, HPX_SIZE_T);
  }
};


/// Object that handles action registration for Node
template <typename Record>
class NodeRegistrar {
 public:
  using node_t = Node<Record>;

  NodeRegistrar() {
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        node_t::partition_node_,
                        node_t::partition_node_handler,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT);
  }
};


/// Object that handles action registration for Tree
template <typename Source, typename Target, typename Record,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class TreeRegistrar {
 public:
  using tree_t = Tree<Source, Target, Record, Expansion, Method>;

  TreeRegistrar() {
    // Tree actions
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::setup_basics_,
                        tree_t::setup_basics_handler,
                        HPX_POINTER, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::delete_tree_,
                        tree_t::delete_tree_handler,
                        HPX_POINTER, HPX_INT, HPX_INT, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        tree_t::recv_node_,
                        tree_t::recv_node_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::send_node_,
                        tree_t::send_node_handler,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::assign_points_,
                        tree_t::assign_points_to_unif_grid,
                        HPX_POINTER, HPX_INT, HPX_POINTER, HPX_INT,
                        HPX_POINTER, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::group_points_,
                        tree_t::group_points_on_unif_grid,
                        HPX_POINTER, HPX_INT, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::merge_points_,
                        tree_t::merge_points_handler,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::merge_points_same_s_and_t_,
                        tree_t::merge_points_same_s_and_t_handler,
                        HPX_POINTER, HPX_INT, HPX_POINTER, HPX_ADDR);
  }
};


/// Object that handles action registration for DualTree
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class DualTreeRegistrar {
 public:
  using dualtree_t = DualTree<Source, Target, Expansion, Method>;

  DualTreeRegistrar() {
    // DualTree actions
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::set_domain_geometry_,
                        dualtree_t::set_domain_geometry_handler,
                        HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        dualtree_t::domain_geometry_init_,
                        dualtree_t::domain_geometry_init_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        dualtree_t::domain_geometry_op_,
                        dualtree_t::domain_geometry_op_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::init_partition_,
                        dualtree_t::init_partition_handler,
                        HPX_ADDR, HPX_ADDR, HPX_INT, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        dualtree_t::recv_points_,
                        dualtree_t::recv_points_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::send_points_,
                        dualtree_t::send_points_handler,
                        HPX_INT, HPX_POINTER, HPX_POINTER, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::create_dual_tree_,
                        dualtree_t::create_dual_tree_handler,
                        HPX_ADDR, HPX_ADDR, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::finalize_partition_,
                        dualtree_t::finalize_partition_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::source_apply_method_,
                        dualtree_t::source_apply_method_handler,
                        HPX_POINTER, HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::source_apply_method_child_done_,
                        dualtree_t::source_apply_method_child_done_handler,
                        HPX_POINTER, HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::target_apply_method_,
                        dualtree_t::target_apply_method_handler,
                        HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::destroy_DAG_LCOs_,
                        dualtree_t::destroy_DAG_LCOs_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::termination_detection_,
                        dualtree_t::termination_detection_handler,
                        HPX_ADDR, HPX_POINTER, HPX_SIZE_T, HPX_POINTER,
                        HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::create_S_expansions_from_DAG_,
                        dualtree_t::create_S_expansions_from_DAG_handler,
                        HPX_ADDR, HPX_POINTER, HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::create_T_expansions_from_DAG_,
                        dualtree_t::create_T_expansions_from_DAG_handler,
                        HPX_ADDR, HPX_POINTER, HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::instigate_dag_eval_,
                        dualtree_t::instigate_dag_eval_handler,
                        HPX_ADDR, HPX_POINTER, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        dualtree_t::instigate_dag_eval_remote_,
                        dualtree_t::instigate_dag_eval_remote_handler,
                        HPX_POINTER, HPX_SIZE_T);
  }
};


} // namespace dashmm


#endif // __DASHMM_REGISTRAR_H__

#ifndef __DASHMM_REGISTRAR_H__
#define __DASHMM_REGISTRAR_H__


#include "dashmm/expansionlco.h"
#include "dashmm/targetlco.h"
#include "dashmm/tree.h"



namespace dashmm {


template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class TargetLCORegistrar {
public:
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method,
                                DistroPolicy>;

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
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        targetlco_t::create_, targetlco_t::create_at_locality,
                        HPX_POINTER, HPX_SIZE_T);
  }
};


template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class ExpansionLCORegistrar {
public:
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method,
                                      DistroPolicy>;

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
                        expansionlco_t::spawn_out_edges_handler,
                        HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::spawn_out_edges_from_remote_,
                        expansionlco_t::spawn_out_edges_from_remote_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::create_from_expansion_,
                        expansionlco_t::create_from_expansion_handler,
                        HPX_POINTER, HPX_SIZE_T);
  }
};


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


template <typename Source, typename Target, typename Record,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class TreeRegistrar {
 public:
  using tree_t = Tree<Source, Target, Record, Expansion, Method, DistroPolicy>;

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


template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class DualTreeRegistrar {
 public:
  using dualtree_t = DualTree<Source, Target, Expansion, Method, DistroPolicy>;

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
                        HPX_POINTER, HPX_SIZE_T);
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
                        dualtree_t::edge_lists_,
                        dualtree_t::edge_lists_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::instigate_dag_eval_,
                        dualtree_t::instigate_dag_eval_handler,
                        HPX_ADDR, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        dualtree_t::instigate_dag_eval_remote_,
                        dualtree_t::instigate_dag_eval_remote_handler,
                        HPX_POINTER, HPX_SIZE_T);
  }
};


} // namespace dashmm


#endif // __DASHMM_REGISTRAR_H__
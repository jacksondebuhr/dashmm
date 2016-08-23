#ifndef __REGISTRAR_H__
#define __REGISTRAR_H__


#include "hpx/hpx.h"

#include "dashmm/arrayref.h"
using dashmm::ArrayRef;

#include "tree.h"


template <typename Record>
class NodeRegistrar {
 public:
  using node_t = Node<Record>;

  NodeRegistrar() {
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        node_t::partition_node_,
                        node_t::partition_node_handler,
                        HPX_POINTER, HPX_POINTER, HPX_INT);
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
                        HPX_POINTER, HPX_INT);
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
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT, HPX_ADDR);
  }
};


template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class Registrar {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion, DistroPolicy>;
  using sourceref_t = ArrayRef<Source>;
  using targetref_t = ArrayRef<Target>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using sourcetree_t = Tree<Source, Target, Source, Expansion, Method,
                            DistroPolicy>;
  using targettree_t = Tree<Source, Target, Target, Expansion, Method,
                            DistroPolicy>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method, DistroPolicy>;
  using distropolicy_t = DistroPolicy;

  Registrar() : sreg_{}, treg_{}, streg_{}, ttreg_{} {
    // DualTree actions
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        dualtree_t::set_domain_geometry_,
                        dualtree_t::set_domain_geometry_handler,
                        HPX_ADDR, HPX_ADDR, HPX_ADDR);
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
                        HPX_ADDR, HPX_ADDR, HPX_INT, HPX_ADDR);
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
  }

private:
  NodeRegistrar<Source> sreg_;
  NodeRegistrar<Target> treg_;
  TreeRegistrar<Source, Target, Source, Expansion, Method, DistroPolicy> streg_;
  TreeRegistrar<Source, Target, Target, Expansion, Method, DistroPolicy> ttreg_;
};

#endif // __REGISTRAR_H__
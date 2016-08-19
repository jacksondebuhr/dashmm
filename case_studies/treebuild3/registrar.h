#ifndef __REGISTRAR_H__
#define __REGISTRAR_H__


#include "hpx/hpx.h"

#include "dashmm/arrayref.h"
using dashmm::ArrayRef;

#include "tree.h"


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
  using sourcenode_t = Node<Source, Target, Source, Expansion, Method,
                            DistroPolicy>;
  using targetnode_t = Node<Source, Target, Target, Expansion, Method,
                            DistroPolicy>;
  using tree_t = DualTree<Source, Target, Expansion, Method, DistroPolicy>;
  using distropolicy_t = DistroPolicy;

  Registrar() {
    // One for each sort of node.
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        sourcenode_t::partition_node_,
                        sourcenode_t::partition_node_handler,
                        HPX_POINTER, HPX_POINTER, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        targetnode_t::partition_node_,
                        targetnode_t::partition_node_handler,
                        HPX_POINTER, HPX_POINTER, HPX_INT);

    // DualTree actions
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::set_domain_geometry_,
                        tree_t::set_domain_geometry_handler,
                        HPX_ADDR, HPX_ADDR, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        tree_t::domain_geometry_init_,
                        tree_t::domain_geometry_init_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        tree_t::domain_geometry_op_,
                        tree_t::domain_geometry_op_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::init_partition_,
                        tree_t::init_partition_handler,
                        HPX_ADDR, HPX_ADDR, HPX_INT, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::merge_points_source_,
                        tree_t::merge_points_handler<Source>,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::merge_points_target_,
                        tree_t::merge_points_handler<Target>,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        tree_t::recv_points_,
                        tree_t::recv_points_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::send_points_,
                        tree_t::send_points_handler,
                        HPX_INT, HPX_POINTER, HPX_POINTER, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        tree_t::recv_node_source_,
                        tree_t::recv_node_handler<Source>,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        tree_t::recv_node_target_,
                        tree_t::recv_node_handler<Target>,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::send_node_source_,
                        tree_t::send_node_handler<Source>,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::send_node_target_,
                        tree_t::send_node_handler<Target>,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::create_dual_tree_,
                        tree_t::create_dual_tree_handler,
                        HPX_ADDR, HPX_ADDR, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::finalize_partition_,
                        tree_t::finalize_partition_handler,
                        HPX_ADDR);
  }
};

#endif // __REGISTRAR_H__
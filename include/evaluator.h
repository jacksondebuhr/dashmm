// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_EVALUATOR_H__
#define __DASHMM_EVALUATOR_H__


#include "include/particle.h"


template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename, typename> class Method>
class Evaluator {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, expansion_t>;

  // Related type aliases to avoid clutter
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using sourcenode_t = SourceNode<Source, Target, Expansion, Method>;
  using targetnode_t = TargetNode<Source, Target, Expansion, Method>;

  Evaluator() {
    // TODO: register all actions in the system

    // TargetLCO related
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        targetlco_t::init_, targetlco_t::init_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        targetlco_t::operation_,
                        targetlco_t::operation_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        targetlco_t::predicate_,
                        targetlco_t::predicate_handler,
                        HPX_POINTER, HPX_SIZE_T);

    // ExpansionLCO related
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        targetlco_t::init_, targetlco_t::init_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        expansionlco_t::operation_,
                        expansionlco_t::operation_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        expansionlco_t::predicate_,
                        expansionlco_t::predicate_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_m_,
                        expansionlco_t::s_to_m_handler,
                        HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE,
                        HPX_DOUBLE, HPX_DOUBLE, HPX_ADDR, HPX_INT, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_l_,
                        expansionlco_t::_s_to_l_handler,
                        HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE,
                        HPX_DOUBLE, HPX_DOUBLE, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::m_to_m_,
                        expansionlco_t::m_to_m_handler,
                        HPX_ADDR, HPX_INT, HPX_INT, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::m_to_l_,
                        expansionlco_t::m_to_l_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::l_to_l_,
                        expansionlco_t::l_to_l_handler,
                        HPX_ADDR, HPX_INT, HPX_INT, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::m_to_t_,
                        expansionlco_t::m_to_t_handler,
                        HPX_INT, HPX_DOUBLE, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::l_to_t_,
                        expansionlco_t::l_to_t_handler,
                        HPX_INT, HPX_DOUBLE, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_t_,
                        expansionlco_t::s_to_t_handler,
                        HPX_POINTER, HPX_INT, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::add_,
                        expansionlco_t::add_handler,
                        HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::create_from_expansion_,
                        expansionlco_t::create_from_expansion_handler,
                        HPX_POINTER, HPX_SIZE_T);

    // Source Node related
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        sourcenode_t::self_delete_,
                        sourcenode_t::self_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        sourcenode_t::node_delete_,
                        sourcenode_t::node_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        sourcenode_t::child_done_,
                        sourcenode_t::child_done_handler,
                        HPX_POINTER, HPX_ADDR, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
                        sourcenode_t::partition_,
                        sourcenode_t::partition_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

    // Target Node related
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        targetnode_t::self_delete_,
                        targetnode_t::self_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        targetnode_t::node_delete_,
                        targetnode_t::node_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
                        targetnode_t::partition_,
                        targetnode_t::partition_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
  }
 private:
};


#endif // __DASHMM_EVALUATOR_H__

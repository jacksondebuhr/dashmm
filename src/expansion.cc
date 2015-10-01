#include "expansion.h"

#include <map>

#include "hpx/hpx.h"
#include "libsync/sync.h"

//other dashmm stuff


namespace dashmm {


//The range of Expansion type identifiers
constexpr int kFirstExpansionType = 1000;
constexpr int kLastExpansionType = 1999;

//A global kept by locality zero
int next_available_expansion_ = kFirstExpansionType;


//map from type to description of that type
std::map<int, ExpansionDesc> expansion_table_;


//parameter types for user-marshalled actions
struct register_expansion_params_t {
  int type;
  ExpansionDesc desc;
  hpx_action_t coregen;
};


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int register_expansion_handler(register_expansion_params_t *parms,
                               size_t size) {
  assert(expansion_table_[parms->type].count() == 0
            && "Registering expansion over existing");
  parms->desc.core_data = nullptr;
  if (parms->desc.core_size) {
    coregen_function_t func = hpx_action_gen_handler(parms->coregen);
    parms->desc.core_data = func(parms->desc.size);
  }
  expansion_table_[parms->type] = parms->desc;
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           register_expansion_action, register_expansion_handler,
           HPX_POINTER, HPX_SIZE_T);


int request_next_expansion_identifier_handler(void *UNUSED) {
  int retval = sync_fadd(&next_available_expansion_, 1, SYNC_ACQ_REL);
  assert(retval <= kLastExpansionType && "Out of expansions, somehow...");
  HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, 0, request_next_expansion_identifier_action,
           request_next_expansion_identifier_handler, HPX_POINTER);


/////////////////////////////////////////////////////////////////////
// Internal Routines
/////////////////////////////////////////////////////////////////////


const ExpansionDesc &get_expansion_desc(int type) {
  auto entry = expansion_table_.find(type);
  if (entry == expansion_table_.end()) {
    return expansion_table_[kFirstExpansionType];
  } else {
    return *entry;
  }
}


int request_next_expansion_identifier() {
  int retval{0};
  hpx_call_sync(HPX_THERE(0), request_next_expansion_identifier_action,
                &retval, sizeof(retval), nullptr);
  return retval;
}


/////////////////////////////////////////////////////////////////////
// Public Interface
/////////////////////////////////////////////////////////////////////


void Expansion::destroy() {
}


Point Expansion::center() const {
}


std::complex<double> Expansion::term(size_t i) const {
}


std::unique_ptr<Expansion> Expansion::S_to_M(Point center,
                                  std::vector<Source>::iterator first,
                                  std::vector<Source>::iterator last) const {
  //
}


std::unique_ptr<Expansion> Expansion::S_to_L(Point center,
                                  std::vector<Source>::iterator first,
                                  std::vector<Source>::iterator last) const {
  //
}


std::unique_ptr<Expansion> Expansion::M_to_M(int from_child,
                                             double s_size) const {
  //
}


std::unique_ptr<Expansion> Expansion::M_to_L(Index s_index, double s_size,
                                  Index t_index) const {
  //
}


std::unique_ptr<Expansion> Expansion::L_to_L(int to_child,
                                             double t_size) const {
}


void Expansion::M_to_T(std::vector<Target>::iterator first,
            std::vector<Target>::iterator last) const {
  //
}


void Expansion::L_to_T(std::vector<Target>::iterator first,
            std::vector<Target>::iterator last) const {
  //
}


void Expansion::S_to_T(std::vector<Source>::iterator s_first,
            std::vector<Source>::iterator s_last,
            std::vector<Target>::iterator t_first,
            std::vector<Target>::iterator t_last) const {
  //
}

void Expansion::add_expansion(const Expansion *temp1) {
}


void Expansion::from_sum(const std::vector<const Expansion *> &exps) {
}


std::unique_ptr<Expansion> Expansion::get_new_expansion(Point center) const {
}


int register_expansion(ExpansionDesc desc, hpx_action_t coregen) {
  register_expansion_params_t parms{};
  parms.type = request_next_expansion_identifier();
  parms.desc = desc;
  parms.coregen = coregen;
  hpx_bcast_rsync(register_expansion_action, &parms, sizeof(parms));
  return parms.type;
}


} // namespace dashmm

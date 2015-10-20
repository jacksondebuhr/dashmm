#include "include/expansionref.h"

#include <cassert>

#include <hpx/hpx.h>

//other dashmm


namespace dashmm {


int ExpansionRef::type() const {
  assert(valid());
  hpx_addr_t from = hpx_addr_add(data_, );
}


bool ExpansionRef::provides_L() const {
  assert(valid());
}


size_t ExpansionRef::size() const {
  assert(valid());
}


Point ExpansionRef::center() const {
  assert(valid());
}


std::complex<double> ExpansionRef::term(size_t i) const {

}


std::unique_ptr<Expansion> ExpansionRef::S_to_M(Point center,
                                std::vector<Source>::iterator first,
                                std::vector<Source>::iterator last) const {
  //
}


std::unique_ptr<Expansion> ExpansionRef::S_to_L(Point center,
                                std::vector<Source>::iterator first,
                                std::vector<Source>::iterator last) const {
  //
}


std::unique_ptr<Expansion> ExpansionRef::M_to_M(int from_child,
                                                double s_size) const {
  //
}


std::unique_ptr<Expansion> ExpansionRef::M_to_L(Index s_index, double s_size,
                                  Index t_index) const {
  //
}


std::unique_ptr<Expansion> ExpansionRef::L_to_L(int to_child,
                                                double t_size) const {
  //
}


void ExpansionRef::M_to_T(std::vector<Target>::iterator first,
            std::vector<Target>::iterator last) const {
  //
}


void ExpansionRef::L_to_T(std::vector<Target>::iterator first,
            std::vector<Target>::iterator last) const {
  //
}


void ExpansionRef::S_to_T(std::vector<Source>::iterator s_first,
            std::vector<Source>::iterator s_last,
            std::vector<Target>::iterator t_first,
            std::vector<Target>::iterator t_last) const {
  //
}


void ExpansionRef::add_expansion(const Expansion *temp1) {

}


void ExpansionRef::from_sum(const std::vector<const Expansion *> &exps) {

}


std::unique_ptr<Expansion> ExpansionRef::get_new_expansion(Point center) const {

}


ExpansionSerial *ExpansionRef::pin() {
  if (data == HPX_NULL) {
    return nullptr;
  }
  ExpansionSerial *retval{nullptr};
  hpx_gas_try_pin(data_, (void **)&retval);
  return retval;
}


void ExpansionRef::unpin() {
  if (data_ == HPX_NULL) {
    return;
  }
  hpx_gas_unpin(data_);
}



ExpansionRef globalize_expansion(Expansion *exp, hpx_addr_t where) {
  if (exp == nullptr) {
    return ExpansionRef{HPX_NULL};
  }
  ExpansionSerialPtr serial = met->serialize();
  size_t size = serial->size + sizeof(ExpansionSerial);
  hpx_addr_t data = hpx_gas_alloc_local_at_sync(size, 0, where);
  assert(data != HPX_NULL);
  hpx_gas_memput_rsync(data, serial.get(), size);
  return data;
}


} // namespace dashmm

#include "include/expansionref.h"

#include <cassert>

#include <hpx/hpx.h>

//other dashmm


namespace dashmm {


//we repeat this so we can use it below.
std::unique_ptr<Expansion> create_expansion(int type, size_t size, void *data)


int ExpansionRef::type() const {
  assert(valid());
  hpx_addr_t from = hpx_addr_add(data_, );
}


bool ExpansionRef::provides_L() const {
  assert(valid());
  ExpansionSerial *local = pin();
  if (local) {
    bool retval = local->provides_L;
    unpin();
    return retval;
  } else {
    assert(0 && "This should not happen. We are in SMP.")
  }
}


size_t ExpansionRef::size() const {
  assert(valid());
  ExpansionSerial *local = pin();
  if (local) {
    size_t retval = local->term_count;
    unpin();
    return retval;
  } else {
    assert(0 && "This should not happen. We are in SMP.")
  }
}


Point ExpansionRef::center() const {
  assert(valid());
  ExpansionSerial *local = pin();
  if (local) {
    Point retval = local->center;
    unpin();
    return retval;
  } else {
    assert(0 && "This should not happen. We are in SMP.")
  }
}


std::complex<double> ExpansionRef::term(size_t i) const {
  setup_local_expansion();
  return exp_->term(i);
}


std::unique_ptr<Expansion> ExpansionRef::S_to_M(Point center,
                                      Source *first, Source *last) const {
  setup_local_expansion();
  return exp_->S_to_M(center, first, last);
}


std::unique_ptr<Expansion> ExpansionRef::S_to_L(Point center,
                                      Source *first, Source *last) const {
  setup_local_expansion();
  return exp_->S_to_L(center, first, last);
}


std::unique_ptr<Expansion> ExpansionRef::M_to_M(int from_child,
                                                double s_size) const {
  setup_local_expansion();
  return exp_->M_to_M(from_child, s_size);
}


std::unique_ptr<Expansion> ExpansionRef::M_to_L(Index s_index, double s_size,
                                  Index t_index) const {
  setup_local_expansion();
  return exp_->M_to_L(s_index, s_size, t_index);
}


std::unique_ptr<Expansion> ExpansionRef::L_to_L(int to_child,
                                                double t_size) const {
  setup_local_expansion();
  return exp_->L_to_L(to_child, t_size);
}


void ExpansionRef::M_to_T(Target *first, Target *last) const {
  setup_local_expansion();
  exp_->M_to_T(first, last);
}


void ExpansionRef::L_to_T(Target *first, Target *last) const {
  setup_local_expansion();
  exp_->L_to_T(first, last);
}


void ExpansionRef::S_to_T(Source *s_first, Source *s_last,
                          Target *t_first, Target *t_last) const {
  setup_local_expansion();
  exp_->S_to_T(s_first, s_last, t_first, t_last);
}


void ExpansionRef::add_expansion(const Expansion *temp1) {
  setup_local_expansion();
  exp_->add_expansion(temp1);
  save_to_global();
}


void ExpansionRef::from_sum(const std::vector<const Expansion *> &exps) {
  setup_local_expansion();
  exp_->from_sum(exps);
  save_to_global();
}


std::unique_ptr<Expansion> ExpansionRef::get_new_expansion(Point center) const {
  setup_local_expansion();
  return exp_->get_new_expansion(center);
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


//NOTE: This makes a copy. Which will need to happen in non SMP mode eventually,
// but for now, this might be a performance hit.
void ExpansionRef::setup_local_expansion() {
  if (exp_ != nullptr) {
    return;
  }
  ExpansionSerial *local = pin();
  size_t total = sizeof(ExpansionSerial) + local->size;
  exp_ = create_expansion(local->type, total, local);
  unpin();
}


//NOTE: There is some repetition between this and globalize_expansion...
void ExpansionRef::save_to_global() {
  assert(exp_ != nullptr);
  ExpansionSerialPtr serial = exp_->serialize();
  size_t size = serial->size + sizeof(ExpansionSerial);
  hpx_gas_memput_rsync(data_, serial.get(), size);
}


ExpansionRef globalize_expansion(Expansion *exp, hpx_addr_t where) {
  if (exp == nullptr) {
    return ExpansionRef{HPX_NULL};
  }
  ExpansionSerialPtr serial = exp->serialize();
  size_t size = serial->size + sizeof(ExpansionSerial);
  hpx_addr_t data = hpx_gas_alloc_local_at_sync(size, 0, where);
  assert(data != HPX_NULL);
  hpx_gas_memput_rsync(data, serial.get(), size);
  return data;
}


} // namespace dashmm

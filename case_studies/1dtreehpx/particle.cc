#include "particle.h"


int Particle_set_approx(Particle *local, double phi, hpx_addr_t whendone) {
  local->set_approx(phi);
  assert(whendone != HPX_NULL);
  hpx_lco_set(whendone, 0, nullptr, HPX_NULL, HPX_NULL);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           Particle_set_approx_action, Particle_set_approx,
           HPX_POINTER, HPX_DOUBLE, HPX_ADDR);


SetApproxContinuation Particle::set_approx_cont() const {
  return SetApproxContinuation(data_, Particle_set_approx_action);
}


bool operator<(const Particle &a, const Particle &b) {
  return a.x() < b.x();
}



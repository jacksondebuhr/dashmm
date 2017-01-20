#ifndef __ONEDTREE_MOMENTS_MIXIN_H__
#define __ONEDTREE_MOMENTS_MIXIN_H__


#include <cstdio>

#include <hpx/hpx.h>


// The identity and operation for moment reductions
extern hpx_action_t moment_reduction_ident;
extern hpx_action_t moment_reduction_op;


// The data for moments
struct Moment {
  double mtot;
  double xcom;
  double Q00;

  Moment() : mtot{0.0}, xcom{0.0}, Q00{0.0} { }
};


// The mixin class
template <typename T>
struct MomentMixIn {
 public:
  MomentMixIn() : M{} { }

  void compile_test() {
    T *backptr = static_cast<T *>(this);
    if (backptr->left && backptr->right && backptr->parts == nullptr) {
      fprintf(stdout, "Hooray!\n");
    }
  }


  void compute_moments(hpx_addr_t sync = HPX_NULL) {
    hpx_addr_t done = sync;
    if (done == HPX_NULL) {
      done = hpx_lco_future_new(sizeof(Moment));
      assert(done != HPX_NULL);
    }

    T *ptr = static_cast<T *>(this);
    hpx_call(HPX_HERE, compute_moments_, done, &ptr);

    if (sync == HPX_NULL) {
      hpx_lco_wait(done);
      hpx_lco_delete_sync(done);
    }
  }


  void compute_local_moments() {
    T *ptr = static_cast<T *>(this);
    auto P = ptr->parts;

    for (int i = 0; i < ptr->count; ++i) {
      M.mtot += P[i].mass;
      M.xcom += P[i].pos * P[i].mass;
    }
    if (M.mtot > 0.0) {
      M.xcom /= M.mtot;
      for (int i = 0; i < ptr->count; ++i) {
        double dx = P[i].pos - M.xcom;
        M.Q00 += 2.0 * P[i].mass * dx * dx;
      }
    }
  }


  static void register_actions() {
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        compute_moments_,
                        compute_moments_handler,
                        HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        save_and_continue_,
                        save_and_continue_handler,
                        HPX_POINTER, HPX_ADDR);
  }


  Moment M;


 private:
  static hpx_action_t compute_moments_;
  static hpx_action_t save_and_continue_;


  static int compute_moments_handler(T * ptr) {
    if (ptr->is_leaf()) {
      ptr->compute_local_moments();
      return hpx_thread_continue(&ptr->M, sizeof(Moment));
    } else {
      int n_in{0};
      if (ptr->left) ++n_in;
      if (ptr->right) ++n_in;

      hpx_addr_t sync = hpx_lco_reduce_new(n_in, sizeof(Moment),
                                           moment_reduction_ident,
                                           moment_reduction_op);

      if (ptr->left) {
        ptr->left->compute_moments(sync);
      }
      if (ptr->right) {
        ptr->right->compute_moments(sync);
      }

      return hpx_call_when_cc(sync, HPX_HERE, save_and_continue_, &ptr, &sync);
    }
  }


  // Action to save the moment and clean up the reducer before continuing
  // the value up the tree.
  static int save_and_continue_handler(T *ptr, hpx_addr_t done) {
    hpx_lco_get(done, sizeof(Moment), &ptr->M);
    hpx_lco_delete_sync(done);
    return hpx_thread_continue(&ptr->M, sizeof(Moment));
  }
};

template <typename T>
hpx_action_t MomentMixIn<T>::compute_moments_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t MomentMixIn<T>::save_and_continue_ = HPX_ACTION_NULL;


#endif // __ONEDTREE_MOMENTS_MIXIN_H__
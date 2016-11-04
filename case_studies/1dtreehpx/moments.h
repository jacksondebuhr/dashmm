#ifndef __ONEDTREE_MOMENTS_MIXIN_H__
#define __ONEDTREE_MOMENTS_MIXIN_H__


#include <cstdio>


struct Moment {
  double mtot;
  double xcom;
  double Q00;

  Moment() : mtot{0.0}, xcom{0.0}, Q00{0.0} { }
};


template <typename T>
struct MomentMixIn {
 public:
  MomentMixIn() : M{} { }

  void compile_test() {
    T *backptr = static_cast<T *>(this);
    if (this->left && this->right && this->particles == nullptr) {
      fprintf(stdout, "Hooray!\n");
    }
  }

  Moment M;
};


#endif // __ONEDTREE_MOMENTS_MIXIN_H__
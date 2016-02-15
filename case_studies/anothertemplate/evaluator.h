#ifndef __EVALUATOR_H__
#define __EVALUATOR_H__


#include "sourceref.h"


template <typename T>
class Evaluator {
 private:
  hpx_action_t the_action_id_;

  static int accumulate_handler(T *vec, size_t bytes) {
    size_t count = bytes / sizeof(T);
    assert(count == 2);
    std::cout << "Accumulate received: " << vec[0] << " and " << vec[1]
              << std::endl;

    SourceRef<T> data{vec[0]};
    data.contribute(vec[1]);
    data.destroy();

    hpx_exit(HPX_SUCCESS);
  }

 public:
  Evaluator() {
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, the_action_id_,
                        accumulate_handler, HPX_POINTER, HPX_SIZE_T);

    // Now register stuff for SourceRef
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        SourceRef<T>::contribute_action_,
                        SourceRef<T>::contribute_handler,
                        HPX_POINTER, HPX_SIZE_T);
  }

  int accumulate(T *vec, size_t num) {
    return hpx_run(&the_action_id_, vec, sizeof(T) * num);
  }
};


#endif // __EVALUATOR_H__

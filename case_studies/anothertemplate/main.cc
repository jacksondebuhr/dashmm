#include <iostream>

#include <hpx/hpx.h>


// For this one we want to work out if there is a way to simultaneously
// register for a related class when the instance is created.
//
// This could either be a nested class (I think this one is easy) or a
// related but separate class.

// In any event, we need to work out how a class that will be used multiple
// times in a dashmm instance will work. The actions cannot be a member of the
// class probably. Or, if they are, they cannot be registered by the class.
// Which means that the DASHMM object will have to register them. So basically
// they are static members of the class. Probably private so the users cannot
// muck with them. And so that means we need to make dashmm or the evaluator
// a friend.


template <typename T>
class Evaluator;


template <typename T>
class SourceRef {
 public:
  SourceRef(T a) {
    assert(hpx_is_active());
    data_ = hpx_gas_alloc_local(sizeof(T), sizeof(T), 0);
    T *local{nullptr};
    assert(hpx_gas_try_pin(data_, (void **)&local));
    *local = a;
    std::cout << "I am a SourceRef created with " << a << std::endl;
    hpx_gas_unpin(data_);
  }

  void destroy() {
    hpx_gas_free_sync(data_);
    data_ = HPX_NULL;
  }

  void contribute(T data) {
    hpx_call_sync(data_, contribute_action_, nullptr, 0, &data, sizeof(data));
  }

 private:
  // NOTE: The action handles must be static here. And thus, they cannot access
  // non-static members of the class. In this case, we get around it by
  // using thread_current_target(), but what if more was needed?
  //
  // Would we really need to pass in all the extra stuff to the action? That
  // would be obnoxious...
  static int contribute_handler(T *data, size_t bytes) {
    assert(bytes == sizeof(T));

    hpx_addr_t target = hpx_thread_current_target();

    T *local{nullptr};
    assert(hpx_gas_try_pin(target, (void **)&local));
    std::cout << "I begin contribute with: " << *local << std::endl;
    *local += *data;
    std::cout << "I end contribute with: " << *local << std::endl;
    hpx_gas_unpin(target);

    return HPX_SUCCESS;
  }

  friend class Evaluator<T>;

  static hpx_action_t contribute_action_;
  hpx_addr_t data_;
};

template <typename T>
hpx_action_t SourceRef<T>::contribute_action_ = HPX_ACTION_NULL;


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



Evaluator<int> dashforint{};
Evaluator<double> dashfordouble{};


int main(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) {
    return -1;
  }

  int argsi[2] = {10, 15};
  dashforint.accumulate(argsi, 2);

   double args[2] = {2.0, 3.14};
   dashfordouble.accumulate(args, 2);

  hpx_finalize();
  return 0;
}

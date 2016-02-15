#ifndef __SOURCE_REF_H__
#define __SOURCE_REF_H__


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

  // NOTE: we declare evaluator a friend so that it might register this class'
  // actions.
  friend class Evaluator<T>;

  static hpx_action_t contribute_action_;
  hpx_addr_t data_;
};

template <typename T>
hpx_action_t SourceRef<T>::contribute_action_ = HPX_ACTION_NULL;


#endif // __SOURCE_REF_H__

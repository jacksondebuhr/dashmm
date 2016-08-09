#ifndef __DASHMM_RANKWISE_H__
#define __DASHMM_RANKWISE_H__


#include <hpx/hpx.h>


extern hpx_action_t rankwise_init_action;


template <typename T>
class RankLocal {
public:
  RankLocal(hpx_addr_t global)
        : local_{nullptr}, global_{global}, remote_{false} {
    if (global_ == HPX_NULL) return;

    if (!hpx_gas_try_pin(global_, (void **)&local_)) {
      //If the pin fails, assume that it is because the address is not local
      // so get that data
      local_ = hpx_malloc_registered(sizeof(T));
      if (local_) {
        remote_ = true;
        hpx_gas_memget_sync(local_, global_, sizeof(T));
      } else {
        global_ = HPX_NULL;
      }
    }
  }

  ~RankLocal() {
    if (local_ != nullptr) {
      if (remote_) {
        hpx_free_registered(local_);
      } else {
        hpx_gas_unpin(global_);
      }
      global_ = HPX_NULL;
      local_ = nullptr;
    }
  }

  RankLocal(const RankLocal &other) = delete;
  RankLocal &operator=(const RankLocal &other) = delete;

  RankLocal(RankLocal &&other) {
    local_ = other.local_;
    global_ = other.global_;
    remote_ = other.remote_;
    other.local_ = nullptr;
    other.global_ = HPX_NULL;
    other.remote_ = false;
  }

  RankLocal &operator=(RankLocal &&other) {
    local_ = other.local_;
    global_ = other.global_;
    remote_ = other.remote_;
    other.local_ = nullptr;
    other.global_ = HPX_NULL;
    other.remote_ = false;
  }

  bool valid() const {return global_ != HPX_NULL && local_ != nullptr;}
  bool remote() const {return remote_;}

  T& operator*() {return *local_;}
  T* operator->() {return local_;}

private:
  T *local_;
  hpx_addr_t global_;
  bool remote_;
};


// OKAY, so this is going to need a constructor from an address
//  The usage pattern seems to be create it on one rank, then share the
//  address via broadcast. At which point, initialization of values is just
//  via normal use of here().
template <typename T>
class RankWise {
 public:
  RankWise(hpx_addr_t dat = HPX_NULL) : data_{dat} { }

  void allocate() {
    assert(data_ == HPX_NULL);
    data_ = hpx_gas_calloc_cyclic(hpx_get_num_ranks(), sizeof(T), 0);
  }

  bool valid() const {return data_ != HPX_NULL;}
  hpx_addr_t data() const {return data_;}

  RankLocal<T> here() const {
    return there(hpx_get_my_rank());
  }

  RankLocal<T> there(int rank) {
    return RankLocal<T>(address_at(rank));
  }

 private:
  hpx_addr_t address_at(int rank) {
    assert(data_ != HPX_NULL);
    return hpx_addr_add(data_, sizeof(T) * rank, sizeof(T));
  }

  hpx_addr_t data_;
};


#endif // __DASHMM_RANKWISE_H__

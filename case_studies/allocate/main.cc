#include <cassert>
#include <cstdio>
#include <cstring>

#include <vector>

#include <hpx/hpx.h>


// The size of the array - this is to mimic the DASHMM use case
constexpr size_t kRecordSize = 48;
constexpr size_t kNumRecords = 2000000;
constexpr size_t kArraySize = kRecordSize * kNumRecords;

// The size of the LCO data
constexpr size_t kKB = 1024;
constexpr size_t kLCOInitSize = kKB;
constexpr size_t kLCOSize = 4 * kKB;
// The number of LCOs - this will scale with the number of particles
constexpr size_t kPrefactor = 10;
constexpr size_t kLCOCount = kPrefactor * kNumRecords / 40;

// The size of the 'local' memory allocation
constexpr size_t kLocalAllocCount = 300000000;

struct Junk {
  char data[24];
};

struct Data {
  int count;
  char payload[];
};


void init_handler(Data *i, size_t bytes, Data *init, size_t init_bytes) {
  i->count = 1;
  memcpy(i->payload, init, init_bytes);
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
           init_action, init_handler,
           HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);


void operation_handler(Data *lhs, const void *rhs, size_t bytes) {
  lhs->count -= 1;
  memcpy(lhs->payload, rhs, bytes);
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
           operation_action, operation_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


bool predicate_handler(const Data *i, size_t bytes) {
  return (i->count == 0);
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
           predicate_action, predicate_handler,
           HPX_POINTER, HPX_SIZE_T);


// We have to forward declare the handler because this is a recursive action
int lco_tree_alloc_handler(size_t idx, hpx_addr_t *lcos, hpx_addr_t done);
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           lco_tree_alloc_action, lco_tree_alloc_handler,
           HPX_SIZE_T, HPX_POINTER, HPX_ADDR);

int lco_tree_set_handler(size_t idx, hpx_addr_t *lcos, hpx_addr_t done);
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           lco_tree_set_action, lco_tree_set_handler,
           HPX_SIZE_T, HPX_POINTER, HPX_ADDR);


// This is a tree spawn that treats the array as if it were an implicit binary
// tree. If needed, we could do this instead as an implicit octtree to better
// match the DASHMM case.
int lco_tree_alloc_handler(size_t idx, hpx_addr_t *lcos, hpx_addr_t done) {
  // create child done detection lco
  hpx_addr_t cdone = hpx_lco_and_new(2);
  assert(cdone != HPX_NULL);

  // spawn work at children if they exist
  size_t l_child = 2 * idx + 1;
  if (l_child < kLCOCount) {
    hpx_call(HPX_HERE, lco_tree_alloc_action, HPX_NULL,
             &l_child, &lcos, &cdone);
  } else {
    hpx_lco_and_set(cdone, HPX_NULL);
  }

  size_t r_child = 2 * idx + 2;
  if (r_child < kLCOCount) {
    hpx_call(HPX_HERE, lco_tree_alloc_action, HPX_NULL,
             &r_child, &lcos, &cdone);
  } else {
    hpx_lco_and_set(cdone, HPX_NULL);
  }

  // allocate my lco
  //lcos[idx] = hpx_lco_future_new(kLCOSize);

  char *initdata = new char[kLCOInitSize];
  lcos[idx] = hpx_lco_user_new(kLCOSize + sizeof(int), init_action,
                               operation_action, predicate_action,
                               initdata, kLCOInitSize);
  assert(lcos[idx] != HPX_NULL);
  delete [] initdata;

  // setup dependent call of done when cdone triggers
  hpx_call_when(cdone, cdone, hpx_lco_delete_action, done, nullptr, 0);

  return HPX_SUCCESS;
}


int lco_tree_set_handler(size_t idx, hpx_addr_t *lcos, hpx_addr_t done) {
  // create child done detection lco
  hpx_addr_t cdone = hpx_lco_and_new(2);
  assert(cdone != HPX_NULL);

  // spawn work at children if they exist
  size_t l_child = 2 * idx + 1;
  if (l_child < kLCOCount) {
    hpx_call(HPX_HERE, lco_tree_set_action, HPX_NULL,
             &l_child, &lcos, &cdone);
  } else {
    hpx_lco_and_set(cdone, HPX_NULL);
  }

  size_t r_child = 2 * idx + 2;
  if (r_child < kLCOCount) {
    hpx_call(HPX_HERE, lco_tree_set_action, HPX_NULL,
             &r_child, &lcos, &cdone);
  } else {
    hpx_lco_and_set(cdone, HPX_NULL);
  }

  // set my lco
  char *data = new char[kLCOSize];
  assert(data != nullptr);
  hpx_lco_set_lsync(lcos[idx], kLCOSize, data, HPX_NULL);
  delete [] data;

  // setup dependent call of done when cdone triggers
  hpx_call_when(cdone, cdone, hpx_lco_delete_action, done, nullptr, 0);

  return HPX_SUCCESS;
}


int rankwise_handler(hpx_addr_t array_alloc_done, hpx_addr_t lco_alloc_done,
                     hpx_addr_t array_clear_done) {
  int rank = hpx_get_my_rank();
  fprintf(stdout, "Hello from rank %d\n", rank);

  // Allocate a large block of memory
  hpx_addr_t array_gas = hpx_gas_alloc_local(1, kArraySize, 0);
  assert(array_gas != HPX_NULL);
  hpx_lco_and_set(array_alloc_done, HPX_NULL);
  fprintf(stdout, "  %d: array allocated\n", rank);
  hpx_lco_wait(array_alloc_done);

  // allocate some local memory
  std::vector<Junk> local_junk{};
  for (size_t i = 0; i < kLocalAllocCount; ++i) {
    //Junk toadd{};
    //toadd.data[14] = 'c';
    local_junk.push_back(Junk{});
  }

  // in parallel allocate a bunch of LCOs
  hpx_addr_t *lcos = new hpx_addr_t[kLCOCount];
  assert(lcos != nullptr);
  hpx_addr_t local_alloc_done = hpx_lco_and_new(1);
  assert(local_alloc_done != HPX_NULL);
  hpx_call_when(local_alloc_done, local_alloc_done, hpx_lco_delete_action,
                lco_alloc_done, nullptr, 0);
  size_t root_idx{0};
  hpx_call(HPX_HERE, lco_tree_alloc_action, HPX_NULL,
           &root_idx, &lcos, &local_alloc_done);
  fprintf(stdout, "  %d: lco allocation spawned\n", rank);
  hpx_lco_wait(lco_alloc_done);

  // set the large chunk to zero
  double *array_local{nullptr};
  assert(hpx_gas_try_pin(array_gas, (void **)&array_local));
  size_t count = kArraySize / sizeof(double);
  for (size_t i = 0; i < count; ++i) {
    array_local[i] = 0.0;
  }
  hpx_lco_and_set(array_clear_done, HPX_NULL);
  fprintf(stdout, "  %d: array cleared\n", rank);
  hpx_lco_wait(array_clear_done);

  // in parallel, start setting the LCOs
  hpx_addr_t local_set_done = hpx_lco_and_new(1);
  hpx_call_when(local_set_done, local_set_done, hpx_lco_delete_action,
                HPX_NULL, nullptr, 0);
  hpx_call(HPX_HERE, lco_tree_set_action, HPX_NULL,
           &root_idx, &lcos, &local_set_done);
  fprintf(stdout, "  %d: array set spawned\n", rank);

  // wait for all of those to finish
  hpx_lco_wait_all(kLCOCount, lcos, nullptr);
  fprintf(stdout, "  %d: array sets done\n", rank);

  // clean up this rank's resources
  hpx_gas_free_sync(array_gas);
  for (size_t i = 0; i < kLCOCount; ++i) {
    hpx_lco_delete_sync(lcos[i]);
  }
  delete [] lcos;

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           rankwise_action, rankwise_handler,
           HPX_ADDR, HPX_ADDR, HPX_ADDR);



int main_handler(int UNUSED) {
  fprintf(stdout, "Hello from main_handler!\n");

  // First allocate the shared synchronization LCOs
  int n_ranks = hpx_get_num_ranks();
  hpx_addr_t array_alloc_done = hpx_lco_and_new(n_ranks);
  hpx_addr_t lco_alloc_done = hpx_lco_and_new(n_ranks);
  hpx_addr_t array_clear_done = hpx_lco_and_new(n_ranks);

  // Then broadcast so each rank will do its work
  hpx_bcast_rsync(rankwise_action, &array_alloc_done, &lco_alloc_done,
                  &array_clear_done);

  // Clean up the shared resources
  hpx_lco_delete_sync(array_alloc_done);
  hpx_lco_delete_sync(lco_alloc_done);
  hpx_lco_delete_sync(array_clear_done);

  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           main_action, main_handler,
           HPX_INT);


int main(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) {
    return -1;
  }

  int unused{7};
  hpx_run(&main_action, nullptr, &unused);

  hpx_finalize();

  return 0;
}

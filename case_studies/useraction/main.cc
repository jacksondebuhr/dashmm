#include <iostream>

#include <hpx/hpx.h>


// TODO: one problem is that the static specifier on all of this really
// fixes this to be a singleton. Not sure how I feel about that. Basically,
// we have to use leaf_ from inside a static member. So perhaps we can make
// all of it non-static?
//
// No we pass leaf in as an argument


template <typename T, typename E>
class ArrayMapAction {
 public:
  typedef void (*leaf_function_t)(T *, const size_t, const E *);

  // Do a default constructor?
  ArrayMapAction(void(*f)(T *, const size_t, const E *)) {
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, leaf_, f,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER);
    if (!registered_) {
      HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, spawn_, spawn_handler,
                          HPX_SIZE_T, HPX_SIZE_T, HPX_ADDR,
                          HPX_ACTION_T, HPX_POINTER, HPX_POINTER);
      HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, root_, root_handler,
                          HPX_SIZE_T, HPX_SIZE_T, HPX_SIZE_T,
                          HPX_ACTION_T, HPX_POINTER, HPX_POINTER);
      registered_ = 1;
    }
  }

  static hpx_action_t root() {return root_;}
  hpx_action_t leaf() const {return leaf_;}
  static int registered() {return registered_;}

 private:
  static int root_handler(size_t count, size_t chunk_count, size_t chunk_size,
                          hpx_action_t leaf, E *env, T *data) {
    hpx_addr_t alldone = hpx_lco_and_new(chunk_count);

    hpx_time_t start = hpx_time_now();

    hpx_call(HPX_HERE, spawn_, HPX_NULL, &count, &chunk_size, &alldone,
             &leaf, &env, &data);

    hpx_lco_wait(alldone);

    hpx_time_t end = hpx_time_now();

    hpx_lco_delete_sync(alldone);

    fprintf(stdout, "Mapping took %lg [us]\n", hpx_time_diff_us(start, end));

    hpx_exit(HPX_SUCCESS);
  }

  static int spawn_handler(size_t count, size_t chunk_size, hpx_addr_t alldone,
                           hpx_action_t leaf, E *env, T *data) {
    if (count <= chunk_size) {
      // fprintf(stdout, "Leaf: %lu %lu\n", count, chunk_size);
      leaf_function_t lfunc = (leaf_function_t)hpx_action_get_handler(leaf);
      lfunc(data, count, env);
      hpx_lco_set_lsync(alldone, 0, nullptr, HPX_NULL);
    } else {
      // split in two and call
      size_t num_chunks = count / chunk_size;
      if (count % chunk_size) {
        ++num_chunks;
      }

      // compute numbers for L and R
      size_t num_left = num_chunks / 2;
      // size_t num_right = num_chunks - num_left;
      size_t count_left = num_left * chunk_size;
      size_t count_right = count - count_left;
      T *data_left = data;
      T *data_right = &data[num_left * chunk_size];

      // fprintf(stdout, "Split: %lu %lu -- %lu %lu\n",
      //        num_left, count_left, num_right, count_right);

      hpx_call(HPX_HERE, spawn_, HPX_NULL, &count_left, &chunk_size, &alldone,
               &leaf, &env, &data_left);
      hpx_call(HPX_HERE, spawn_, HPX_NULL, &count_right, &chunk_size, &alldone,
               &leaf, &env, &data_right);
    }

    return HPX_SUCCESS;
  }

  static int registered_;
  static hpx_action_t root_;
  static hpx_action_t spawn_;
  hpx_action_t leaf_;
};

template <typename T, typename E>
int ArrayMapAction<T, E>::registered_ = 0;

template <typename T, typename E>
hpx_action_t ArrayMapAction<T, E>::root_ = HPX_ACTION_NULL;

template <typename T, typename E>
hpx_action_t ArrayMapAction<T, E>::spawn_ = HPX_ACTION_NULL;


template <typename T, typename E>
int mapper(const ArrayMapAction<T, E> &act, size_t total, size_t over_factor,
           E *env, T *data) {
  size_t n_per_chunk = total / (hpx_get_num_threads() * over_factor);
  size_t n_chunks{0};

  if (n_per_chunk == 0) {
    n_per_chunk = 1;
    n_chunks = total;
  } else {
    // Guarantees that all chunks are the same size. There could be a
    // chunk with less than n_per_chunk at the end.
    n_chunks = total / n_per_chunk;
    if (total % n_per_chunk) {
      ++n_chunks;
    }
  }

  // fprintf(stdout,
  //        "\nIn mapper: count = %lu - n_chunks = %lu - n_per_chunk = %lu\n",
  //        total, n_chunks, n_per_chunk);

  hpx_action_t rt = act.root();
  hpx_action_t leaf = act.leaf();
  return hpx_run(&rt, &total, &n_chunks, &n_per_chunk, &leaf, &env, &data);
}


// Specific examples
struct SourceData {
  double blah;
};

struct Environment {
  double yoink;
};

// an example function
void mappedaction_handler(SourceData *source, const size_t count,
                          const Environment *env) {
  for (size_t i = 0; i < count; ++i) {
    source[i].blah += env->yoink;
  }
}

void boogieboogieboogie(SourceData *source, const size_t count,
                        const Environment *env);


// The user must create an instance of the object, so that it might be passed
// to the array object.
ArrayMapAction<SourceData, Environment> yoink_update(mappedaction_handler);
ArrayMapAction<SourceData, Environment> boogie_update(boogieboogieboogie);


constexpr int kArraySize = 40;


int main(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) {
    return -1;
  }

  // Get some data
  SourceData array[kArraySize];
  for (int i = 0; i < kArraySize; ++i) {
    array[i].blah = i;
  }

  /*
     for (int i = 0; i < kArraySize; ++i) {
     fprintf(stdout, "%lg ", array[i].blah);
     }
     fprintf(stdout, "\n");
     //*/


  // Call the mapper
  Environment env{13.2};
  mapper(yoink_update, kArraySize, 2, &env, array);
  mapper(boogie_update, kArraySize, 2, &env, array);

  /*
     for (int i = 0; i < kArraySize; ++i) {
     fprintf(stdout, "%lg ", array[i].blah);
     }
     fprintf(stdout, "\n");
     //*/


  hpx_finalize();
  return 0;
}


void boogieboogieboogie(SourceData *source, const size_t count,
                        const Environment *env) {
  for (size_t i = 0; i < count; ++i) {
    fprintf(stdout, "boogie: %lg\n", source[i].blah);
  }
}

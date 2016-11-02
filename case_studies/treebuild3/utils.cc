#include "utils.h"


// NOTE: No cleanup needed, this is just a helper for the test driver.


int point_count(char scaling, int n) {
  if (scaling == 'w') {
    return n;
  } else {
    int my_rank = hpx_get_my_rank();
    int num_ranks = hpx_get_num_ranks();
    int nper = n / num_ranks;
    int remain = n % num_ranks;
    return (my_rank < remain ? nper + 1 : nper);
  }
}


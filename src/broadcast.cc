#include "dashmm/broadcast.h"


namespace dashmm {


int broadcast_value_handler(char *value, size_t size) {
  hpx_exit(size, value);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           broadcast_value_action, broadcast_value_handler,
           HPX_POINTER, HPX_SIZE_T);


} // dashmm

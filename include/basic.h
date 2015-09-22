#ifndef __DASHMM_BASIC_INTERFACE_H__
#define __DASHMM_BASIC_INTERFACE_H__


#include "include/builtins.h"
#include "include/types.h"


namespace dashmm {


int init(int *argc, char ***argv);


int finalize();


int evaluate(ObjectHandle source, size_t source_count, ...);
int evaluate(ObjectHandle source, ...);


int allocate_array(size_t count, size_t size, ObjectHandle *obj);


int deallocate_array(ObjectHandle obj);


int array_put(ObjectHandle obj, size_t first, size_t last,
                     size_t size, void *in_data);


int array_get(ObjectHandle obj, size_t first, size_t last,
                     size_t size, void *out_data);


} // namespace dashmm


#endif // __DASHMM_BASIC_INTERFACE_H__

#ifndef __DASHMM_DAG_TO_JSON_H__
#define __DASHMM_DAG_TO_JSON_H__


#include <vector>

#include "dashmm/daginfo.h"


namespace dashmm {


void output_dag_as_JSON(std::vector<DAGNode *> &source,
                        std::vector<DAGNode *> &target,
                        std::vector<DAGNode *> &internal);


} // dashmm


#endif

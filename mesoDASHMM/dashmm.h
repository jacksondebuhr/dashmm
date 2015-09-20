#ifndef __DASHMM_H__
#define __DASHMM_H__

#include "method.h"


namespace dashmm {


Method *bh_method(double crit);

Expansion *bh_expansion();

Method *fmm_method();

Expansion *fmm_expansion(const int p);


void init();

void evaluate(std::vector<Source> &sources,
              std::vector<Target> &targets,
              double limit,
              Method *method,
              Expansion *expand);

void finalize();


} //namespace dashmm


#endif // __DASHMM_H__

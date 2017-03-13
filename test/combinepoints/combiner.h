#ifndef __DASHMM_TEST_COMBINEPOINTS_COMBINER_H__
#define __DASHMM_TEST_COMBINEPOINTS_COMBINER_H__


#include <string>

#include "../common/common.h"


class Combiner {
public:
  Combiner() : fname_{}, head_{}, sources_{nullptr}, targets_{nullptr} { }

  // Can throw std::ios_base::failure
  Combiner(std::string fname);

  ~Combiner() {
    if (sources_) {
      delete [] sources_;
    }
    if (targets_) {
      delete [] targets_;
    }
  }

  // This will throw std::runtime_error
  void Combine(const Combiner &other);

  // This will throw std::ios_base::failure
  void Write(std::string fname);

  std::string FileName() const {return fname_;}

  bool HasAllPhis() const {
    return head_.has_laplace && head_.has_yukawa && head_.has_helmholtz;
  }

private:
  std::string fname_;
  FileHeader head_;
  FileSourceData *sources_;
  FileTargetData *targets_;
}



#endif // __DASHMM_TEST_COMBINEPOINTS_COMBINER_H__
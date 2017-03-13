#include "combiner.h"

#include <cstdio>

#include <functional>
#include <ios>
#include <stdexcept>

Combiner::Combiner(std::string fname)
    : fname_{fname},
      head_{},
      sources_{nullptr},
      targets_{nullptr} {
  FILE *ofd = fopen(fname_.c_str(), "rb");
  if (ofd == nullptr) {
    std::string message{"Error opening file: "};
    message.append(fname_);
    throw std::ios_base::failure(message);
  }

  if (1 != fread(&head_, sizeof(FileHeader), 1, ofd)) {
    fclose(ofd);
    std::string message{"Error reading header in file: "};
    message.append(fname_);
    throw std::ios_base::failure(message);
  }

  if (head_.n_sources) {
    sources_ = new FileSourceData[head_.n_sources];
    if ((unsigned int)head_.n_sources
        != fread(sources_, sizeof(FileSourceData), head_.n_sources, ofd)) {
      delete [] sources_;
      sources_ = nullptr;
      fclose(ofd);
      std::string message{"Error reading sources in file: "};
      message.append(fname_);
      throw std::ios_base::failure(message);
    }
  }

  if (head_.n_targets) {
    targets_ = new FileTargetData[head_.n_targets];
    if ((unsigned int)head_.n_targets
        != fread(targets_, sizeof(FileTargetData), head_.n_targets, ofd)) {
      delete [] sources_;
      sources_ = nullptr;
      delete [] targets_;
      targets_ = nullptr;
      fclose(ofd);
      std::string message{"Error reading targets in file: "};
      message.append(fname_);
      throw std::ios_base::failure(message);
    }
  }

  fclose(ofd);
}

void Combiner::Combine(const Combiner &other) {
  // That's quite the type...
  std::function<void(FileTargetData &a, const FileTargetData &b)>
      copy_relevant_field{nullptr};

  // do some simple checks, and some simple setup
  if (head_.n_sources != other.head_.n_sources ||
      head_.n_targets != other.head_.n_targets) {
    std::string message{"Source and target mismatch between '"};
    message.append(this->FileName());
    message.append("' and '");
    message.append(other.FileName());
    message.append("'.");
    throw std::runtime_error(message);
  }
  if (other.head_.has_laplace) {
    if (head_.has_laplace) {
      std::string message{"Files both contain laplace data: '"};
      message.append(this->FileName());
      message.append("' and '");
      message.append(other.FileName());
      message.append("'.");
      throw std::runtime_error(message);
    } else {
      head_.has_laplace = true;
      copy_relevant_field = [](FileTargetData &a, const FileTargetData &b) {
        a.phi_laplace = b.phi_laplace;
      };
    }
  } else if (other.head_.has_yukawa) {
    if (head_.has_yukawa) {
      std::string message{"Files both contain yukawa data: '"};
      message.append(this->FileName());
      message.append("' and '");
      message.append(other.FileName());
      message.append("'.");
      throw std::runtime_error(message);
    } else {
      head_.has_yukawa = true;
      copy_relevant_field = [](FileTargetData &a, const FileTargetData &b) {
        a.phi_yukawa = b.phi_yukawa;
      };
    }
  } else if (other.head_.has_helmholtz) {
    if (head_.has_helmholtz) {
      std::string message{"Files both contain helmholtz data: '"};
      message.append(this->FileName());
      message.append("' and '");
      message.append(other.FileName());
      message.append("'.");
      throw std::runtime_error(message);
    } else {
      head_.has_helmholtz = true;
      copy_relevant_field = [](FileTargetData &a, const FileTargetData &b) {
        a.phi_helmholtz = b.phi_helmholtz;
      };
    }
  }

  // loop over each target in this
  //  add the data from the equivalent target in other (are they the same
  //  index? They should be.)
  for (int i = 0; i < head_.n_targets; ++i) {
    if (targets_[i].index == other.targets_[i].index) {
      copy_relevant_field(targets_[i], other.targets_[i]);
    } else {
      std::string message{"Index mismatch in files: '"};
      message.append(this->FileName());
      message.append("' and '");
      message.append(other.FileName());
      message.append("'.");
      throw std::runtime_error(message);
    }
  }
}

void Combiner::Write(std::string fname) {
  FILE *ofd = fopen(fname.c_str(), "wb");
  if (ofd == nullptr) {
    std::string message{"Error opening file: "};
    message.append(fname);
    throw std::ios_base::failure(message);
  }

  if (1 != fwrite(&head_, sizeof(head_), 1, ofd)) {
    fclose(ofd);
    std::string message{"Error writing header to file: "};
    message.append(fname_);
    throw std::ios_base::failure(message);
  }

  if ((unsigned int)head_.n_sources
        != fwrite(sources_, sizeof(FileSourceData), head_.n_sources, ofd)) {
    fclose(ofd);
    std::string message{"Error writing sources to file: "};
    message.append(fname_);
    throw std::ios_base::failure(message);
  }

  if ((unsigned int)head_.n_targets
        != fwrite(targets_, sizeof(FileTargetData), head_.n_targets, ofd)) {
    fclose(ofd);
    std::string message{"Error writing targets to file: "};
    message.append(fname_);
    throw std::ios_base::failure(message);
  }
}
// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#include "user_expansion.h"

#include <cstdio>


//We import only what we are using from the dashmm namespace. Details of
// Point and Target are in comments below. Expansion is just the abstract
// base from which User is derived. For details on Point and Index, please
// see ../../include/point.h and ../../include/index.h.
using dashmm::Expansion;
using dashmm::Point;
using dashmm::Source;
using dashmm::Index;
using dashmm::Target;
using dashmm::dcomplex_t;


// This is the specific detail of the implementation of the data for the
// User expansion. There are some requirements on the serialization of an
// expansion, and they are represented below. The start of the serialization
// must be three integers. The first is reserved for use by DASHMM. The second
// is the type identifier of the Expansion. The third is the accuracy parameter
// for this instance of the Expansion. After these, the data can take any form
// so long as it is contiguous.
struct UserData {
  int reserved;
  int type;
  int acc;
  //Whatever else that is needed. One example might be
  //
  // int count;
  // double terms[];
  //
  // For a real number expansion with count terms. The form of this
  // extra data will be dictated by the details of the expansion.
  //
  // Another example is that one might require an offset into a table of
  // precomputed values. (For more on this, see register_user_with_dashmm()
  // below.)
  //
  // For now, we settle for the center point.
  Point center;
};


// Currently, DASHMM expects that the data provided by release() was allocated
// using malloc().
User::User(Point center, int n_digits) {
  // If there was more complication to UserData, this next line would need to
  // be modified.
  bytes_ = sizeof(UserData);
  acc_ = n_digits;
  data_ = static_cast<UserData *>(malloc(bytes_));
  assert(valid());
  data_->type = type();
  data_->acc = acc_;
  data_->center = center;
}


// If the object is valid, free the data when this object is destroyed.
User::~User() {
  if (valid()) {
    free(data_);
    data_ = nullptr;
  }
}


// Note that not only does release provide the data to the caller, but it will
// render this object invalid. Otherwise, when this object was destroyed, an
// attempt would be made to free the data, which is precisely what is not
// intended.
void *User::release() {
  UserData *retval = data_;
  data_ = nullptr;
  return retval;
}


// Currently, User does not have any terms. Typically, the result of this
// function will depend on acc_.
size_t User::size() const {
  return 0;
}

Point User::center() const {
  if (valid()) {
    return data_->center;
  } else {
    return Point(0.0, 0.0, 0.0);
  }
}


// User currently has no terms, so we return a default constructed complex
// number. This would be modified for any real use-case.
dcomplex_t User::term(size_t i) const {return dcomplex_t{};}


// The following operations (S->M, S->L, M->M, M->L, L->L, M->T, L->T, S->T)
// would all have to be implemented for the the specific use-case. For now,
// these all have a simple implementation that spits out a message for each
// operation.
//
// Source is a structure that holds the relevant source point information. The
// relevant portions of its declaration are reproduced here.
//
// struct Source {
//   double charge;
//   Point position;
// };
//
// Similarly, Target holds target data:
//
// struct Target {
//   Point position;
//   dcomplet_t phi;
// }
//

void User::S_to_M(Point center, Source *first, Source *last,
                  double scale) const {
  fprintf(stdout, "S->M for %ld sources\n", last - first);
}

std::unique_ptr<Expansion> User::S_to_L(Point center,
                                          Source *first, Source *last,
                                          double scale) const {
  fprintf(stdout, "S->L for %ld sources\n", last - first);
  return std::unique_ptr<Expansion>{new User{center, acc_}};
}

std::unique_ptr<Expansion> User::M_to_M(int from_child,
                                          double s_size) const {
  fprintf(stdout, "M->M from child %d\n", from_child);
  double h = s_size / 2;
  double px = data_->center.x() + (from_child % 2 == 0 ? h : -h);
  double py = data_->center.y() + (from_child % 4 <= 1 ? h : -h);
  double pz = data_->center.z() + (from_child < 4 ? h : -h);
  return std::unique_ptr<Expansion>{new User{Point{px, py, pz}, acc_}};
}

std::unique_ptr<Expansion> User::M_to_L(Index s_index, double s_size,
                                           Index t_index) const {
  fprintf(stdout, "M->L\n");
  int t2s_x = s_index.x() - t_index.x();
  int t2s_y = s_index.y() - t_index.y();
  int t2s_z = s_index.z() - t_index.z();
  double tx = data_->center.x() - t2s_x * s_size;
  double ty = data_->center.y() - t2s_y * s_size;
  double tz = data_->center.z() - t2s_z * s_size;
  return std::unique_ptr<Expansion>{new User{Point{tx, ty, tz}, acc_}};
}

std::unique_ptr<Expansion> User::L_to_L(int to_child,
                                          double t_size) const {
  fprintf(stdout, "L->L to child %d\n", to_child);
  double h = t_size / 2;
  double cx = data_->center.x() + (to_child % 2 == 0 ? -h : h);
  double cy = data_->center.y() + (to_child % 4 <= 1 ? -h : h);
  double cz = data_->center.z() + (to_child < 4 ? -h : h);
  return std::unique_ptr<Expansion>{new User{Point{cx, cy, cz}, acc_}};
}

void User::M_to_T(Target *first, Target *last, double scale) const {
  fprintf(stdout, "M->T for %ld targets\n", last - first);
}

void User::L_to_T(Target *first, Target *last, double scale) const {
  fprintf(stdout, "L->T for %ld targets\n", last - first);
}

void User::S_to_T(Source *s_first, Source *s_last,
                  Target *t_first, Target *t_last) const {
  fprintf(stdout, "S->T for %ld sources and %ld targets\n",
          s_last - s_first, t_last - t_first);
}

void User::add_expansion(const Expansion *temp1) {
  fprintf(stdout, "Adding an expansion\n");
}


// This will create a new object with this as a prototype.
//
// Note that DASHMM expects these objects to be created with new.
std::unique_ptr<Expansion> User::get_new_expansion(Point center) const {
  return std::unique_ptr<Expansion>{new User{center, acc_}};
}


// This is the implementation of the Expansion creation function. Again,
// DASHMM expects the object to be created with new.
//
// Note the macro after the function definition. This registers the function
// with the HPX-5 runtime system. The first argument indicates that this is
// registered as an HPX_FUNCTION. The second argument is 0 because there are no
// relevant options for HPX_FUNCTIONs. Next comes an identifier for this
// function. This identifier will need to be used in
// dashmm::register_user_expansion(). Next is the name of the function being
// registered. Finally are a list of argument types to the function.
//
// When implementing your Expansion, the only things that need to be changed
// in that macro use would be the name of the function, and the name of the
// identifier. It is useful that the function and the identifier are related,
// but that is not required.
Expansion *create_user_expansion(double cx, double cy, double cz,
                                         int acc) {
  Expansion *retval = new User{Point{cx, cy, cx}, acc};
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           create_user_expansion_action, create_user_expansion,
           HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE, HPX_INT);


// This is the implementation of the Expansion interpretation function.
// Note again that this object is created with new.
//
// Also note the HPX-5 registration macro. The overall use is quite similar.
// Once again, the only thing that needs to be changed in this is the
// function name and the chosen identifier. This identifier is also used
// when the Expansion is registered via dashmm::register_user_expansion().
Expansion *interpret_user_expansion(void *data, size_t bytes, int acc) {
  Expansion *retval = new User{(UserData *)data, bytes, acc};
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           interpret_user_expansion_action, interpret_user_expansion,
           HPX_POINTER, HPX_SIZE_T, HPX_INT);


// If there is a table of values required by the expansion, one good place to
// initialize that table is here. Of course, other options are possible, but
// this is convenient. It could be that this function would have to be extended
// to accept an accuracy parameter, but that is a simple change.
void register_user_with_dashmm() {
  dashmm::ReturnCode ret = dashmm::register_user_expansion(
      kUserExpansionType, create_user_expansion_action,
      interpret_user_expansion_action);
  assert(ret == dashmm::kSuccess);
}

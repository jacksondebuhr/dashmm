#include "include/laplace_com.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>


namespace dashmm {


LaplaceCOM::LaplaceCOM(Point center) {
  bytes_ = sizeof(LaplaceCOMData);
  data_ = malloc(bytes_);
  assert(valid());
  data_->type = type();
}

LaplaceCOM::~LaplaceCOM() {
  if (valid()) {
    free(data_);
    data_ = nullptr;
  }
}


std::complex<double> LaplaceCOM::term(size_t i) const {
  if (i == 0) {
    return std::complex<double>{data_->mtot};
  } else if (i < 4) {
    return std::complex<double>{data_->xcom[i - 1]};
  } else {
    return std::complex<double>{data_->Q[i - 4]};
  }
}


std::unique_ptr<Expansion> LaplaceCOM::S_to_M(Point center,
                                              Source *first,
                                              Source *last) const {
  LaplaceCOM *retval{new LaplaceCOM{center}};
  assert(retval);
  retval->calc_mtot(first, last);
  retval->calc_xcom(first, last);
  retval->calc_Q(first, last);
  return std::unique_ptr<Expansion>{retval};
}


std::unique_ptr<Expansion> LaplaceCOM::M_to_M(int from_child,
                                              double s_size) const {
  assert(valid());
  LaplaceCOM *temp = new LaplaceCOM(Point{0.0, 0.0, 0.0});
  temp->set_mtot(data_->mtot);
  temp->set_xcom(data_->xcom);
  temp->set_Q(data_->Q);
  return std::unique_ptr<Expansion>{temp};
}


void LaplaceCOM::M_to_T(Target *first, Target *last) const {
  assert(valid());
  for (auto i = first; i != last; ++i) {
    Point pos{i->position()};

    double sum = 0;
    double diff[3]{pos.x() - data_->xcom[0], pos.y() - data_->xcom[1],
                   pos.z() - data_->xcom[2]};
    double diff2mag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};
    double diffmag{sqrt(diff2mag)};
    double quaddenom{-1.0 / (2.0 * diff2mag * diff2mag * diffmag)};

    sum = data_->mtot / diffmag;
    double qsum{data_->Q[0] * diff[0] * diff[0]};
    qsum += 2.0 * data_->Q[1] * diff[0] * diff[1];
    qsum += 2.0 * data_->Q[2] * diff[0] * diff[2];
    qsum += data_->Q[3] * diff[1] * diff[1];
    qsum += 2.0 * data_->Q[4] * diff[1] * diff[2];
    qsum += data_->Q[5] * diff[2] * diff[2];
    qsum *= quaddenom;
    sum += qsum;

    i->set_phi(i->phi() + std::complex<double>{sum});
  }
}


void LaplaceCOM::S_to_T(Source *s_first, Source *s_last,
                        Target *t_first, Target *t_last) const {
  for (auto targ = t_first; targ != t_last; ++targ) {
    Point pos = targ->position();
    double sum{0.0};
    for (auto i = s_first; i != s_last; ++i) {
      double diff[3]{pos.x() - i->x(), pos.y() - i->y(), pos.z() - i->z()};
      double mag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};
      if (mag > 0) {
        sum += i->charge() / sqrt(mag);
      }
    }

    targ->set_phi(targ->phi() + std::complex<double>{sum});
  }
}


void LaplaceCOM::add_expansion(const Expansion *temp1) {
  double M2 = temp1->term(0).real();
  double D2[3] = {temp1->term(1).real(), temp1->term(2).real()
                  temp1->term(3).real()};
  double Q2[6] = {temp1->term(4).real(), temp1->term(5).real()
                  temp1->term(6).real(), temp1->term(7).real()
                  temp1->term(8).real(), temp1->term(9).real()};

  double Mprime = data->mtot + M2;

  double Dprime[3]{};
  Dprime[0] = (data->mtot * data->xcom[0] + M2 * D2[0]) / Mprime;
  Dprime[1] = (data->mtot * data->xcom[1] + M2 * D2[1]) / Mprime;
  Dprime[2] = (data->mtot * data->xcom[2] + M2 * D2[2]) / Mprime;

  double diff1[3] = {data->xcom[0] - Dprime[0], data->xcom[1] - Dprime[1],
                     data->xcom[2] - Dprime[2]};
  double diff1mag = diff1[0] * diff1[0] + diff1[1] * diff1[1]
                    + diff1[2] * diff1[2];
  double diff2[3] = {D2[0] - Dprime[0], D2[1] - Dprime[1], D2[2] - Dprime[2]};
  double diff2mag = diff2[0] * diff2[0] + diff2[1] * diff2[1]
                    + diff2[2] * diff2[2];

  double Qprime[6]{};
  Qprime[0] = data->mtot * (3.0 * diff1[0] * diff1[0] - diff1mag) + data->Q[0];
  Qprime[1] = data->mtot * (3.0 * diff1[0] * diff1[1]) + data->Q[1];
  Qprime[2] = data->mtot * (3.0 * diff1[0] * diff1[2]) + data->Q[2];
  Qprime[3] = data->mtot * (3.0 * diff1[1] * diff1[1] - diff1mag) + data->Q[3];
  Qprime[4] = data->mtot * (3.0 * diff1[1] * diff1[2]) + data->Q[4];
  Qprime[5] = data->mtot * (3.0 * diff1[2] * diff1[2] - diff1mag) + data->Q[5];
  Qprime[0] += M2 * (3.0 * diff2[0] * diff2[0] - diff2mag) +Q2[0];
  Qprime[1] += M2 * (3.0 * diff2[0] * diff2[1]) + Q2[1];
  Qprime[2] += M2 * (3.0 * diff2[0] * diff2[2]) + Q2[2];
  Qprime[3] += M2 * (3.0 * diff2[1] * diff2[1] - diff2mag) + Q2[3];
  Qprime[4] += M2 * (3.0 * diff2[1] * diff2[2]) + Q2[4];
  Qprime[5] += M2 * (3.0 * diff2[2] * diff2[2] - diff2mag) + Q2[5];

  set_mtot(Mprime);
  set_xcom(Dprime);
  set_Q(Qprime);
}


void LaplaceCOM::calc_mtot(Source *first, Source *last) {
  assert(valid());
  data_->mtot_ = 0.0;
  for (auto i = first; i != last; ++i) {
    data_->mtot_ += i->charge();
  }
}


void LaplaceCOM::calc_xcom(Source *first, Source *last) {
  assert(valid());
  data_->xcom_[0] = 0.0;
  data_->xcom_[1] = 0.0;
  data_->xcom_[2] = 0.0;
  for (auto i = first; i != last; ++i) {
    data_->xcom_[0] += i->charge() * i->x();
    data_->xcom_[1] += i->charge() * i->y();
    data_->xcom_[2] += i->charge() * i->z();
  }
  if (data_->mtot_ != 0.0) {
    double oomtot = 1.0 / mtot_;
    data_->xcom_[0] *= oomtot;
    data_->xcom_[1] *= oomtot;
    data_->xcom_[2] *= oomtot;
  }
}


void LaplaceCOM::calc_Q(Source *first, Source *last) {
  assert(valid());
  for (int i = 0; i < 6; ++i) {
    data_->Q_[i] = 0;
  }

  for (auto i = first; i != last; ++i) {
    double diff[3]{i->x() - data_->xcom_[0], i->y() - data_->xcom_[1],
                   i->z() - data_->xcom_[2]};
    double diffmag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};

    data_->Q_[0] += i->charge() * (3.0 * diff[0] * diff[0] - diffmag);   //Qxx
    data_->Q_[1] += i->charge() * 3.0 * diff[0] * diff[1];               //Qxy
    data_->Q_[2] += i->charge() * 3.0 * diff[0] * diff[2];               //Qxz
    data_->Q_[3] += i->charge() * (3.0 * diff[1] * diff[1] - diffmag);   //Qyy
    data_->Q_[4] += i->charge() * 3.0 * diff[1] * diff[2];               //Qyz
    data_->Q_[5] += i->charge() * (3.0 * diff[2] * diff[2] - diffmag);   //Qzz
  }
}


} //namespace dashmm

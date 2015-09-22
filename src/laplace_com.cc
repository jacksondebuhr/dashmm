#include "include/laplace_com.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>


namespace dashmm {


struct LaplaceCOMData {
  double mtot_;
  double xcom_[3];
  //Ordering: Qxx Qxy Qxz Qyy Qyz Qzz
  double Q_[6];
};


LaplaceCOM::LaplaceCOM(Point center) {
  //(c)allocate space in GAS for it
  //save address
}


void LaplaceCOM::destroy() {
  //free the global memory
  data_ = HPX_NULL;
}


std::complex<double> LaplaceCOM::term(size_t i) const {
  //TODO: if valid, pull the particular value and return it

  if (i == 0) {
    return std::complex<double>{mtot_};
  } else if (i < 4) {
    return std::complex<double>{xcom_[i - 1]};
  } else {
    return std::complex<double>{Q_[i - 4]};
  }
}


std::unique_ptr<Expansion> LaplaceCOM::S_to_M(Point center,
                       std::vector<Source>::iterator first,
                       std::vector<Source>::iterator last) const {
  LaplaceCOM *retval{new LaplaceCOM{center}};
  assert(retval);
  retval->calc_mtot(first, last);
  retval->calc_xcom(first, last);
  retval->calc_Q(first, last);
  return std::unique_ptr<Expansion>{retval};
}


std::unique_ptr<Expansion> LaplaceCOM::M_to_M(int from_child,
                                                double s_size) const {
  LaplaceCOM *temp = new LaplaceCOM(Point{0.0, 0.0, 0.0});
  temp->set_mtot(mtot_);
  temp->set_xcom(xcom_);
  temp->set_Q(Q_);
  return std::unique_ptr<Expansion>{temp};
}


void LaplaceCOM::M_to_T(std::vector<Target>::iterator first,
                          std::vector<Target>::iterator last) const {
  for (auto i = first; i != last; ++i) {
    Point pos{i->position()};

    double sum = 0;
    double diff[3]{pos.x() - xcom_[0], pos.y() - xcom_[1], pos.z() - xcom_[2]};
    double diff2mag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};
    double diffmag{sqrt(diff2mag)};
    double quaddenom{-1.0 / (2.0 * diff2mag * diff2mag * diffmag)};

    sum = mtot_ / diffmag;
    double qsum{Q_[0] * diff[0] * diff[0]};
    qsum += 2.0 * Q_[1] * diff[0] * diff[1];
    qsum += 2.0 * Q_[2] * diff[0] * diff[2];
    qsum += Q_[3] * diff[1] * diff[1];
    qsum += 2.0 * Q_[4] * diff[1] * diff[2];
    qsum += Q_[5] * diff[2] * diff[2];
    qsum *= quaddenom;
    sum += qsum;

    i->set_phi(i->phi() + std::complex<double>{sum});
  }
}


void LaplaceCOM::S_to_T(std::vector<Source>::iterator s_first,
                          std::vector<Source>::iterator s_last,
                          std::vector<Target>::iterator t_first,
                          std::vector<Target>::iterator t_last) const {
  for (auto targ = t_first; targ != t_last; ++targ) {
    Point pos = targ->position();
    double sum{0.0};
    for (auto i = s_first; i != s_last; ++i) {
      double diff[3]{pos.x() - i->x(), pos.y() - i->y(), pos.z() - i->z()};
      double mag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};
      if (mag > 0)
        sum += i->charge() / sqrt(mag);
    }

    targ->set_phi(targ->phi() + std::complex<double>{sum});
  }
}


void LaplaceCOM::add_expansion(const Expansion *temp1) {
  //This version does nothing for LaplaceCOM
}


void LaplaceCOM::from_sum(const std::vector<const Expansion *> &exps) {
  //NOTE: we could add some checks that make sure these expansions are at
  // least all the right size, but we will assume that nobody is up to any
  // shenanigans.
  mtot_ = 0.0;
  xcom_[0] = 0.0;
  xcom_[1] = 0.0;
  xcom_[2] = 0.0;
  for (auto i = exps.begin(); i != exps.end(); ++i) {
    double m_i = (*i)->term(0).real();
    mtot_ += m_i;
    xcom_[0] += m_i * (*i)->term(1).real();
    xcom_[1] += m_i * (*i)->term(2).real();
    xcom_[2] += m_i * (*i)->term(3).real();
  }
  if (mtot_ != 0.0) {
    xcom_[0] /= mtot_;
    xcom_[1] /= mtot_;
    xcom_[2] /= mtot_;
  }

  Q_[0] = 0.0;
  Q_[1] = 0.0;
  Q_[2] = 0.0;
  Q_[3] = 0.0;
  Q_[4] = 0.0;
  Q_[5] = 0.0;
  for (auto i = exps.begin(); i != exps.end(); ++i) {
    double m_i = (*i)->term(0).real();
    double xcom_i[3]{(*i)->term(1).real(),
                     (*i)->term(2).real(), (*i)->term(3).real()};
    double Q_i[6]{(*i)->term(4).real(), (*i)->term(5).real(),
                  (*i)->term(6).real(), (*i)->term(7).real(),
                  (*i)->term(8).real(), (*i)->term(9).real()};

    double diff[3]{xcom_i[0] - xcom_[0],
                   xcom_i[1] - xcom_[1], xcom_i[2] - xcom_[2]};
    double mag2{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};

    Q_[0] += m_i * (3.0 * diff[0] * diff[0] - mag2) + Q_i[0];
    Q_[1] += m_i * 3.0 * diff[0] * diff[1] + Q_i[1];
    Q_[2] += m_i * 3.0 * diff[0] * diff[2] + Q_i[2];
    Q_[3] += m_i * (3.0 * diff[1] * diff[1] - mag2) + Q_i[3];
    Q_[4] += m_i * 3.0 * diff[1] * diff[2] + Q_i[4];
    Q_[5] += m_i * (3.0 * diff[2] * diff[2] - mag2) + Q_i[5];
  }
}



void LaplaceCOM::calc_mtot(std::vector<Source>::iterator first,
                             std::vector<Source>::iterator last) {
  mtot_ = 0.0;
  for (auto i = first; i != last; ++i) {
    mtot_ += i->charge();
  }
}


void LaplaceCOM::calc_xcom(std::vector<Source>::iterator first,
                             std::vector<Source>::iterator last) {
  xcom_[0] = 0.0;
  xcom_[1] = 0.0;
  xcom_[2] = 0.0;
  for (auto i = first; i != last; ++i) {
    xcom_[0] += i->charge() * i->x();
    xcom_[1] += i->charge() * i->y();
    xcom_[2] += i->charge() * i->z();
  }
  if (mtot_ != 0.0) {
    double oomtot = 1.0 / mtot_;
    xcom_[0] *= oomtot;
    xcom_[1] *= oomtot;
    xcom_[2] *= oomtot;
  }
}


void LaplaceCOM::calc_Q(std::vector<Source>::iterator first,
                          std::vector<Source>::iterator last) {
  for (int i = 0; i < 6; ++i) {
    Q_[i] = 0;
  }

  for (auto i = first; i != last; ++i) {
    double diff[3]{i->x() - xcom_[0], i->y() - xcom_[1], i->z() - xcom_[2]};
    double diffmag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};

    Q_[0] += i->charge() * (3.0 * diff[0] * diff[0] - diffmag);   //Qxx
    Q_[1] += i->charge() * 3.0 * diff[0] * diff[1];               //Qxy
    Q_[2] += i->charge() * 3.0 * diff[0] * diff[2];               //Qxz
    Q_[3] += i->charge() * (3.0 * diff[1] * diff[1] - diffmag);   //Qyy
    Q_[4] += i->charge() * 3.0 * diff[1] * diff[2];               //Qyz
    Q_[5] += i->charge() * (3.0 * diff[2] * diff[2] - diffmag);   //Qzz
  }
}

} //namespace dashmm

// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


/// \file src/special_function.cc
/// \brief Implementation of special functions

#include <algorithm>
#include "builtins/special_function.h"


namespace dashmm {


const double M_SQRTPI = 0.9189385332046727417803297L;


void legendre_Plm(int n, double x, double *P) {
  double u = -sqrt(1.0 - x * x);
  P[midx(0, 0)] = 1.0;
  for (int i = 1; i <= n; i++) {
    P[midx(i, i)] = P[midx(i - 1, i - 1)] * u * (2 * i - 1);
  }

  for (int i = 0; i < n; i++) {
    P[midx(i + 1, i)] = P[midx(i, i)] * x * (2 * i + 1);
  }

  for (int m = 0; m <= n; m++) {
    for (int ell = m + 2; ell <= n; ell++) {
      P[midx(ell, m)] = ((2.0 * ell - 1) * x * P[midx(ell - 1, m)] -
                         (ell + m - 1) * P[midx(ell - 2, m)]) / (ell - m);
    }
  }
}

void legendre_Plm_gt1_scaled(int nb, double x, double scale, double *P) {
  double v = scale * x;
  double w = scale * scale;
  double u = sqrt(x * x - 1.0) * scale;

  P[midx(0, 0)] = 1.0;
  for (int n = 1; n <= nb; ++n) {
    P[midx(n, n)] = P[midx(n - 1, n - 1)] * u * (2 * n - 1);
  }

  for (int n = 0; n < nb; ++n) {
    P[midx(n + 1, n)] = P[midx(n, n)] * (2 * n + 1) * v;
  }

  for (int m = 0; m <= nb; ++m) {
    for (int n = m + 2; n <= nb; ++n) {
      P[midx(n, m)] = ((2.0 * n - 1) * v * P[midx(n - 1, m)]
                       - (n + m - 1) * w * P[midx(n - 2, m)]) / (n - m);
    }
  }
}

double Gamma(double x) {
  // the largest argument for which gamma(x) is representable
  const double xbig = 171.624;

  // the smallest positive floating-point number such that 1/xminin is
  // machine representable
  const double xminin = 2.23e-308;

  // the smallest positive floating-point number such that 1 + eps > 1
  const double eps = 2.22e-16;

  // the largest machine representable floating-point number
  const double xinf = 1.79e308;

  // numerator and denominator coefficients for rational minimax
  // approximation over (1, 2)
  const double p[8] = {-1.71618513886549492533811e+0,
                       2.47656508055759199108314e+1,
                       -3.79804256470945635097577e+2,
                       6.29331155312818442661052e+2,
                       8.66966202790413211295064e+2,
                       -3.14512729688483675254357e+4,
                       -3.61444134186911729807069e+4,
                       6.64561438202405440627855e+4};

  const double q[8] = {-3.08402300119738975254353e+1,
                       3.15350626979604161529144e+2,
                       -1.01515636749021914166146e+3,
                       -3.10777167157231109440444e+3,
                       2.25381184209801510330112e+4,
                       4.75584627752788110767815e+3,
                       -1.34659959864969306392456e+5,
                       -1.15132259675553483497211e+5};

  // coefficients for minimax approximation over (12, inf)
  const double c[7] = {-1.910444077728e-03,
                       8.4171387781295e-04,
                       -5.952379913043012e-04,
                       7.93650793500350248e-04,
                       -2.777777777777681622553e-03,
                       8.333333333333333331554247e-02,
                       5.7083835261e-03};

  bool parity = false;
  double y, fact, res;

  fact = 1;
  y = x;

  if (y <= 0.0) {
    // argument is negative
    y = -x;
    double y1 = floor(y);
    res = y - y1;
    if (res != 0.0) {
      if (y1 != floor(y1 * 0.5) * 2) {
        parity = true;
      }
      fact = -M_PI / sin(M_PI * res);
      y += 1.0;
    } else {
      return xinf;
    }
  }

  // argument is positive
  if (y < eps) {
    if (y >= xminin) {
      res = 1.0 / y;
    } else {
      return xinf;
    }
  } else if (y < 12.0) {
    double y1 = y;
    double z;
    int n = 0;
    if (y < 1.0) {
      z = y;
      y += 1.0;
    } else {
      n = floor(y) - 1;
      y -= n;
      z = y - 1.0;
    }

    double xnum = 0.0;
    double xden = 1.0;

    for (int i = 0; i < 8; ++i) {
      xnum = (xnum + p[i]) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;

    if (y1 < y) {
      res = res / y1;
    } else if (y1 > y) {
      for (int i = 0; i < n; ++i) {
        res *= y;
        y += 1.0;
      }
    }
  } else {
    if (y <= xbig) {
      double ysq = y * y;
      double sum = c[6];
      for (int i = 0; i < 6; ++i) {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + M_SQRTPI;
      sum += (y - 0.5) * log(y);
      res = exp(sum);
    } else {
      return xinf;
    }
  }

  if (parity) {
    res = -res;
  }
  if (fact != 1.0) {
    res = fact / res;
  }
  return res;
}

int bessel_In(int nb, double alpha, double x, int ize, double *B) {
  // decimal significance desired
  const int nsig = 16;

  // upper limit on the magnitude of x when ize is 1
  const double exparg = 709.0;

  // upper limit on the magnitude of x when ize is 2
  const double xlarge = 1e4;

  // 10^k, where k is the largest integer such that enten is
  // machine-representable in working precision
  const double enten = 1e308;

  // 10^nsig
  const double ensig = 1e16;

  // 10^(-k) for the smallest integer k such that k >= nsig / 4
  const double rtnsig = 1e-4;

  // smallest abs(x) such that x / 4 does not underflow
  const double enmten = 8.9e-308;

  // Return if nb, alpha, x, or ize are out of range
  if (nb <= 0 || x < 0 || alpha < 0 || alpha >= 1 ||
      (ize == 1 && x > exparg) ||
      (ize == 2 && x > xlarge)) {
    return -1;
  }

  int ncalc = nb;
  double tempa, tempb, tempc;
  double empal, emp2al, halfx, tover, test;
  double p, plast, pold, psave, psavel, sum, en, em;
  int magx, nbmx, n, nstart, nend;

  if (x >= rtnsig) {
    magx = (int) x;
    nbmx = nb - magx;
    n = magx + 1;
    en = (double)(n + n) + (alpha + alpha);
    plast = 1.0;
    p = en / x;

    // Calculate general significance test
    test = ensig + ensig;
    test = (2 * magx > 5 * nsig ? sqrt(test * p) :
            test / pow(1.585, magx));

    bool stepin = true;

    if (nbmx >= 3) {
      // Calculate p-sequence until n = nb - 1. Check for possible overflow
      tover = enten / ensig;
      nstart = magx + 2;
      nend = nb - 1;
      for (int k = nstart; k <= nend; ++k) {
        n = k;
        en += 2.0;
        pold = plast;
        plast = p;
        p = en * plast / k + pold;
        if (p > tover) {
          // To avoid overflow, divide p-sequence by tover. Calculate p-sequence
          // until abs(p) > 1
          tover = enten;
          p /= tover;
          plast /= tover;
          psave = p;
          psavel = plast;
          nstart = n + 1;
          do {
            n += 1;
            en += 2.0;
            pold = plast;
            plast = p;
            p = en * plast / x + pold;
          } while (p <= 1.0);
          tempb = en / x;

          // Calculate backward test, and find ncalc, the highest n such that
          // the test is passed
          test = pold * plast / ensig;
          test = test * (0.5 - 0.5 / (tempb * tempb));
          p = plast * tover;
          n -= 1;
          en -= 2.0;
          nend = (nb <= n ? nb : n);
          for (int ell = nstart; ell <= nend; ++ell) {
            ncalc = ell;
            pold = psavel;
            psavel = psave;
            psave = en * psavel / x + pold;
            if (psave * psavel > test)
              break;
          }
          if (ncalc != nend) {
            ncalc--;
          }

          stepin = false;
          break;
        }
      }

      if (stepin) {
        n = nend;
        en = (double)(n + n) + (alpha + alpha);
        // Calculate special significance test for nbmx >= 3
        test = fmax(test, sqrt(plast * ensig) * sqrt(p + p));
      }
    }

    if (stepin) {
      // Calculate p-sequence
      do {
        n += 1;
        en += 2.0;
        pold = plast;
        plast = p;
        p = en * plast / x + pold;
      } while (p < test);
    }

    // Initialize the backward recursion and the normalization sum
    n += 1;
    en += 2.0;
    tempb = 0.0;
    tempa = 1.0 / p;
    em = (double) n - 1.0;
    empal = em + alpha;
    emp2al = (em - 1.0) + (alpha + alpha);
    sum = tempa * empal * emp2al / em;
    nend = n - nb;

    if (nend < 0) {
      // Store B[n - 1] and set higher orders to 0.0
      B[n - 1] = tempa;
      nend = -nend;
      for (int ell = 0; ell < nend; ++ell) {
        B[n + ell] = 0.0;
      }
    } else {
      // Recur backward via difference equation, calculating (but not storing
      // b[n - 1]), until n = nb
      for (int ell = 0; ell < nend; ++ell) {
        n -= 1;
        en -= 2.0;
        tempc = tempb;
        tempb = tempa;
        tempa = en * tempb / x + tempc;
        em -= 1.0;
        emp2al -= 1.0;
        if (n == 1) {
          break;
        }
        if (n == 2) {
          emp2al = 1.0;
        }
        empal -= 1.0;
        sum = (sum + tempa * empal) * emp2al / em;
      }

      // Store B[nb - 1]
      B[n - 1] = tempa;
      if (nb <= 1) {
        sum = (sum + sum) + tempa;
        if (alpha != 0.0) {  // statement 230
          sum = sum * Gamma(alpha + 1.0) * pow(x * 0.5, -alpha);
        }
        if (ize == 1) {
          sum *= exp(-x);
        }
        tempa = enmten;
        if (sum > 1.0) {
          tempa *= sum;
        }

        for (int n = 0; n < nb; ++n) {
          if (B[n] < tempa) {
            B[n] = 0.0;
          }
          B[n] /= sum;
        }

        return ncalc;
      }

      // Calculate and store B[nb - 2]
      n -= 1;
      en -= 2.0;
      B[n - 1] = (en * tempa) / x + tempb;
      if (n == 1) {
        sum = (sum + sum) + B[0];
        if (alpha != 0.0) { // statement 230
          sum = sum * Gamma(alpha + 1.0) * pow(x * 0.5, -alpha);
        }
        if (ize == 1) {
          sum *= exp(-x);
        }
        tempa = enmten;
        if (sum > 1.0) {
          tempa *= sum;
        }

        for (int n = 0; n < nb; ++n) {
          if (B[n] < tempa) {
            B[n] = 0.0;
          }
          B[n] /= sum;
        }

        return ncalc;
      }
      em -= 1.0;
      emp2al -= 1.0;
      if (n == 2) {
        emp2al = 1.0;
      }
      empal -= 1.0;
      sum = (sum + B[n - 1] * empal) * emp2al / em;
    }

    nend = n - 2;

    if (nend > 0) {
      // Calculate via difference equation and store B[n - 1], until n = 2
      for (int ell = 0; ell < nend; ++ell) {
        n -= 1;
        en -= 2.0;
        B[n - 1] = (en * B[n]) / x + B[n + 1];
        em -= 1.0;
        emp2al -= 1.0;
        if (n == 2) {
          emp2al = 1.0;
        }
        empal -= 1.0;
        sum = (sum + B[n - 1] * empal) * emp2al / em;
      }
    }

    // Calculate B[0]
    B[0] = 2.0 * empal * B[1] / x + B[2];
    sum = (sum + sum) + B[0];

    // Normalize, divide all B[n] by sum
    if (alpha != 0.0) {  // statement 230
      sum = sum * Gamma(alpha + 1.0) * pow(x * 0.5, -alpha);
    }
    if (ize == 1) {
      sum *= exp(-x);
    }
    tempa = enmten;
    if (sum > 1.0) {
      tempa *= sum;
    }

    for (int n = 0; n < nb; ++n) {
      if (B[n] < tempa) {
        B[n] = 0.0;
      }
      B[n] /= sum;
    }
  } else {
    empal = 1.0 + alpha;
    halfx = (x > enmten ? x * 0.5 : 0.0);
    tempa = (alpha != 0 ? pow(halfx, alpha) / Gamma(empal) : 1.0);
    if (ize == 2) {
      tempa *= exp(-x);
    }
    tempb = ((x + 1.0) > 1.0 ? halfx * halfx : 0.0);
    B[0] = tempa + tempa * tempb / empal;
    if (x != 0 && B[0] == 0.0) {
      ncalc = 0;
    }

    if (nb > 1) {
      if (x == 0.0) {
        for (int i = 1; i < nb; ++i) {
          B[i] = 0.0;
        }
      } else {
        // Calculate higher order functions
        tempc = halfx;
        tover = (tempb != 0 ? enmten / tempb :
                 (enmten + enmten) / x);

        for (int n = 1; n < nb; ++n) {
          tempa /= empal;
          empal += 1.0;
          tempa *= tempc;
          if (tempa <= tover * empal) {
            tempa = 0.0;
          }
          B[n] = tempa + tempa * tempb / empal;
          if (B[n] == 0 && ncalc > n) {
            ncalc = n - 1;
          }
        }
      }
    }
  }

  return ncalc;
}

void bessel_in_scaled(int nb, double x, double scale, double *B) {
  const double ensig = 1e-4;
  const double enmten = 1e-300;

  if (x <= ensig) {
    // Use 2-term Taylor expansion
    double scale_x = x / scale;
    double t1 = 1.0;
    double t2 = 0.5 * x * x;
    B[0] = t1 * (1.0 + t2 / 3.0);
    for (int i = 1; i <= nb; ++i) {
      t1 = t1 * scale_x / (2 * i + 1);
      if (t1 <= enmten) {
        t1 = 0.0;
      }
      B[i] = t1 * (1.0 + t2 / (2 * i + 3));
    }
  } else if (x > 1.0e2) {
    for (int  i = 0; i <= nb; ++i) {
      B[i] = 0.0;
    }
  } else {
    double factor = sqrt(M_PI_2 / x);

    assert(bessel_In(nb + 1, 0.5, x, 1, B) == (nb + 1));

    for (int i = 0; i <= nb; ++i) {
      B[i] *= factor;
      factor /= scale;
      if (fabs(B[i]) <= enmten) {
        factor = 0.0;
      }
    }
  }
}

void bessel_kn_scaled(int nb, double x, double scale, double *B) {
  const double xmin = 4.46e-308;
  const double xinf = 1.79e308;
  const double xlarge = 1e8;
  int ncalc = 0;

  if (nb >= 0 && x >= xmin && x < xlarge) {
    double ex = x;
    double p = exp(-ex) / ex * M_PI_2;
    B[0] = p;
    B[1] = p * scale * (1.0 + 1.0 / ex);
    ncalc = 1;
    double u1 = scale / ex;
    double u2 = scale * scale;
    for (int n = 2; n <= nb; ++n) {
      if (fabs(B[n - 1]) * u1 >= xinf / (2 * n - 1)) {
        break;
      }
      B[n] = (2 * n - 1) * u1 * B[n - 1] + u2 * B[n - 2];
      ncalc++;
    }

    for (int n = ncalc + 1; n <= nb; ++n) {
      B[n] = 0.0;
    }
  } else {
    B[0] = 0.0;
    ncalc = (nb + 1 <= 0 ? nb + 1 : 0) - 1;
  }

  assert(ncalc == nb);
}


} // namespace dashmm


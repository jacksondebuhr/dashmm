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

void legendre_Plm_evan_scaled(int nb, double x, double scale, double *P) {
  double v = scale * x;
  double w = scale * scale;
  double u = sqrt(w + v * v);

  P[midx(0, 0)] = 1;

  for (int n = 1; n <= nb; ++n) {
    P[midx(n, n)] = P[midx(n - 1, n - 1)] * u * (2 * n - 1);
  }

  for (int n = 0; n < nb; ++n) {
    P[midx(n + 1, n)] = P[midx(n, n)] * (2 * n + 1) * v;
  }

  for (int m = 0; m <= nb; ++m) {
    for (int n = m + 2; n <= nb; ++n) {
      P[midx(n, m)] = ((2.0 * n - 1) * v * P[midx(n - 1, m)]
                       + (n + m - 1) * w * P[midx(n - 2, m)]) / (n - m);
    }
  }
}

void legendre_Plm_prop_scaled(int nb, double x, double scale, dcomplex_t *P) {
  double u = sqrt(1 - x * x);
  dcomplex_t v{0.0, scale};
  double s2 = scale * scale;

  P[midx(0, 0)] = 1;

  for (int n = 1; n <= nb; ++n)
    P[midx(n, n)] = P[midx(n - 1, n - 1)] * v * u * (2.0 * n - 1);

  for (int n = 0; n < nb; ++n)
    P[midx(n + 1, n)] = (2.0 * n + 1) * x * (-v) * P[midx(n, n)];

  for (int m = 0; m <= nb; ++m) {
    for (int n = m + 2; n <= nb; ++n) {
      P[midx(n, m)] =
        (-v * (2.0 * n - 1) * x * P[midx(n - 1, m)]
         + s2 * (n + m - 1.0) * P[midx(n - 2, m)]) / (n * 1.0 - m);
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

int bessel_Jn(int nb, double alpha, double x, double *B) {
  const double ENTEN = 1e308;
  const double ENSIG = 1e17;
  const double RTNSIG = 1e-4;
  const double ENMTEN = 8.9e-308;
  const double XLARGE = 1e4;
  const double TWOPI1 = 6.28125;
  const double TWOPI2 = 1.935307179586476925286767e-3;
  const double PI2 = 0.636619772367581343075535;
  const double FACT[25] =
    {1.0, 1.0, 2.0, 6.0, 24.0, 1.2e2, 7.2e2, 5.04e3,
     4.032e4, 3.6288e5, 3.6288e6, 3.99168e7, 4.790016e8, 6.2270208e9,
     8.71782912e10, 1.307674368e12, 2.0922789888e13, 3.55687428096e14,
     6.402373705728e15, 1.21645100408832e17, 2.43290200817664e18,
     5.109094217170944e19, 1.12400072777760768e21,
     2.585201673888497664e22, 6.2044840173323943936e23};

  // Return if nb, alpha, x are out of range
  if (nb <= 0 || x < 0 || x > XLARGE || alpha < 0 || alpha >= 1) {
    return -1;
  }

  int I, J, K, L, M, MAGX, N, NBMX, NCALC, NEND, NSTART;
  double ALPEM, ALP2EM, CAPP, CAPQ, EM, EN, GNU, HALFX, P,
    PLAST, POLD, PSAVE, PSAVEL, S, SUM, T, T1, TEMPA, TEMPB, TEMPC,
    TEST, TOVER, XC, XIN, XK, XM, VCOS, VSIN, Z;

  MAGX = (int) x;
  NCALC = nb;
  for (I = 0; I < nb; ++I)
    B[I] = 0.0;

  if (x < RTNSIG) {
    // Two-term ascending series for small x.
    TEMPA = 1.0;
    ALPEM = 1.0 + alpha;
    HALFX = 0.0;
    if (x > ENMTEN)
      HALFX = 0.5 * x;
    if (alpha != 0.0)
      TEMPA = pow(HALFX, alpha) / (alpha * Gamma(alpha));
    TEMPB = 0.0;
    if ((x + 1.0) > 1.0)
      TEMPB = -HALFX * HALFX;
    B[0] = TEMPA + TEMPA * TEMPB / ALPEM;
    if (x != 0.0 && B[0] == 0.0)
      NCALC = 0;
    if (nb != 1) {
      if (x <= 0.0) {
        for (N = 1; N < nb; ++N)
          B[N] = 0.0;
      } else {
        // Calculate higher order functions
        TEMPC = HALFX;
        TOVER = (ENMTEN + ENMTEN) / x;
        if (TEMPB != 0.0)
          TOVER = ENMTEN / TEMPB;
        for (N = 2; N <= nb; ++N) {
          TEMPA /= ALPEM;
          ALPEM += 1.0;
          TEMPA *= TEMPC;
          if (TEMPA <= TOVER * ALPEM)
            TEMPA = 0.0;
          B[N - 1] = TEMPA + TEMPA * TEMPB / ALPEM;
          if (B[N - 1] == 0.0 && NCALC > N)
            NCALC = N - 1;
        }
      }
    }
  } else if (x > 25.0 && nb <= MAGX + 1) {
    // Asymptotic series for x > 25.0
    XC = sqrt(PI2 / x);
    XIN = pow(0.125 / x, 2);
    M = 11;
    if (x >= 35.0)
      M = 8;
    if (x >= 130.0)
      M = 4;
    XM = 4.0 * M;
    // Argument reduction for sin and cos routines
    T = (int) (x / (TWOPI1 + TWOPI2) + 0.5);
    Z = ((x - T * TWOPI1) - T * TWOPI2) - (alpha + 0.5) / PI2;
    VSIN = sin(Z);
    VCOS = cos(Z);
    GNU = alpha + alpha;

    for (I = 1; I <= 2; ++I) {
      S = ((XM - 1.0) - GNU)  * ((XM - 1.0) + GNU) * XIN * 0.5;
      T = (GNU - (XM - 3.0)) * (GNU + (XM - 3.0));
      CAPP = S * T / FACT[2 * M];
      T1 = (GNU - (XM + 1.0)) * (GNU + (XM + 1.0));
      CAPQ = S * T1 / FACT[2 * M + 1];
      XK = XM;
      K = M + M;
      T1 = T;
      for (J = 2; J <= M; ++J) {
        XK -= 4.0;
        S = ((XK - 1.0) - GNU) * ((XK - 1.0) + GNU);
        T = (GNU - (XK - 3.0)) * (GNU + (XK - 3.0));
        CAPP = (CAPP + 1.0 / FACT[K - 2]) * S * T * XIN;
        CAPQ = (CAPQ + 1.0 / FACT[K - 1]) * S * T1 * XIN;
        K -= 2;
        T1 = T;
      }
      CAPP += 1.0;
      CAPQ = (CAPQ + 1.0) * (GNU * GNU - 1.0) * (0.125 / x);
      B[I - 1] = XC * (CAPP * VCOS - CAPQ * VSIN);
      if (nb == 1) { // goto 300
        return NCALC;
      }

      T = VSIN;
      VSIN = -VCOS;
      VCOS = T;
      GNU += 2.0;
    }

    // If nb > 2, compute j(x, order +i), i = 2, nb - 1
    if (nb > 2) {
      GNU = alpha + alpha + 2.0;
      for (J = 3; J <= nb; ++J) {
        B[J - 1] = GNU * B[J - 2] / x - B[J - 3];
        GNU += 2.0;
      }
    }
  } else {
    // User recurrence to generate results. First initialize the calculation of
    // P * S.
    NBMX = nb - MAGX;
    N = MAGX + 1;
    EN = (double) (N + N) + (alpha + alpha);
    PLAST = 1.0;
    P = EN / x;

    // Calculate general significance test
    TEST = ENSIG + ENSIG;
    bool stepin = true;

    if (NBMX >= 3) {
      // Calculate P * S until N = nb - 1. Check for possible overflow.
      TOVER = ENTEN / ENSIG;
      NSTART = MAGX + 2;
      NEND = nb - 1;
      EN = (double) (NSTART + NSTART) - 2.0 + (alpha + alpha);
      for (K = NSTART; K <= NEND; ++K) {
        N = K;
        EN += 2.0;
        POLD = PLAST;
        PLAST = P;
        P = EN * PLAST / x - POLD;
        if (P > TOVER) {
          // To avoid overflow, divie P * S by TOVER. Calculate P * S until
          // abs(P) > 1
          TOVER = ENTEN;
          P /= TOVER;
          PLAST /= TOVER;
          PSAVE = P;
          PSAVEL = PLAST;
          NSTART = N + 1;
          do {
            N += 1;
            EN += 2.0;
            POLD = PLAST;
            PLAST = P;
            P = EN * PLAST / x - POLD;
          } while (P <= 1.0);
          TEMPB = EN / x;

          // Calculate backward test and find NCALC, the highest n such that the
          // test is passed.
          TEST = POLD * PLAST * (0.5 - 0.5 / (TEMPB * TEMPB));
          TEST = TEST / ENSIG;
          P = PLAST * TOVER;
          N -= 1;
          EN -= 2.0;
          NEND = (nb <= N ? nb : N);

          for (L = NSTART; L <= NEND; ++L) {
            POLD = PSAVEL;
            PSAVEL = PSAVE;
            PSAVE = EN * PSAVEL / x - POLD;
            if (PSAVE * PSAVEL > TEST) {
              break;
              // NCALC = L - 1; (FORTRAN)
              // goto 190 (FORTRAN)
            }
          }
          //NCALC = NEND; (FORTRAN)

          // In the original Fortran code (1) The if-branch is entered, NCALC is
          // set to L - 1 and then jumps to statement 190. (2) The if-branch is
          // not entered, L is NEND + 1 exiting the for-loop, and NCALC is set
          // to NEND, which is L - 1. So we set NCALC to L - 1 in both cases.
          NCALC = L - 1;

          // And flip boolean stepin to false and exits the for-loop on K.
          stepin = false;
          break; // goto 190
        } // P > TOVER
      } // K = NSTART:NEND

      if (stepin) {
        N = NEND;
        EN = (double) (N + N) + (alpha + alpha);
        // Calculate special significance test for NBMX > 2
        TEST = std::max(TEST, sqrt(PLAST * ENSIG) * sqrt(P + P));
      }
    } // NBMX >= 3

    if (stepin) {
      // Calculate P * S until significance test passes.
      do {
        N += 1;
        EN += 2.0;
        POLD = PLAST;
        PLAST = P;
        P = EN * PLAST / x - POLD;
      } while (P < TEST);
    }

    // Initialize the backward recursion and the normalization sum
    N += 1;     // statement 190
    EN += 2.0;
    TEMPB = 0.0;
    TEMPA = 1.0 / P;
    M = 2 * N - 4 * (N / 2);
    SUM = 0;
    EM = (double) (N / 2);
    ALPEM = (EM - 1.0) + alpha;
    ALP2EM = (EM + EM) + alpha;
    if (M != 0)
      SUM = TEMPA * ALPEM * ALP2EM / EM;
    NEND = N - nb;

    if (NEND > 0) {
      // Recur backward via difference equation, calculating (but not storing)
      // B[N - 1], until N = nb.
      for (L = 1; L <= NEND; ++L) {
        N -= 1;
        EN -= 2.0;
        TEMPC = TEMPB;
        TEMPB = TEMPA;
        TEMPA = (EN * TEMPB) / x - TEMPC;
        M = 2 - M;
        if (M != 0) {
          EM -= 1.0;
          ALP2EM = (EM + EM) + alpha;
          if (N == 1) {
            // goto 210 (FORTRAN)
            break;
          }
          ALPEM = (EM - 1.0) + alpha;
          if (ALPEM == 0.0)
            ALPEM = 1.0;
          SUM = (SUM + TEMPA * ALP2EM) * ALPEM / EM;
        }
      }
    } // NEND > 0

    // Store B[nb - 1]
    B[N - 1] = TEMPA;  // statement 210
    if (NEND >= 0) {
      if (nb <= 1) {
        ALP2EM = alpha;
        if ((alpha + 1) == 1.0)
          ALP2EM = 1.0;
        SUM = SUM + B[0] * ALP2EM;
        // goto 250
        // Normalize. Divide all B[N] by SUM.
        if ((alpha + 1.0) != 1.0)
          SUM = SUM * Gamma(alpha) * pow(x * 0.5, -alpha);
        TEMPA = ENMTEN;
        if (SUM > 1.0)
          TEMPA *= SUM;
        for (N = 0; N < nb; ++N) {
          if (fabs(B[N]) < TEMPA)
            B[N] = 0.0;
          B[N] = B[N] / SUM;
        }
        return NCALC;
      } else {
        // Calculate and store B[nb - 2]
        N -= 1;
        EN -= 2.0;
        B[N - 1] = (EN * TEMPA) / x - TEMPB;
        if (N == 1) {
          //goto 240
          EM -= 1; // statement 240
          ALP2EM = (EM + EM) + alpha;
          if (ALP2EM == 0)
            ALP2EM = 1;
          SUM += B[0] * ALP2EM;

          // Normalize. Divide all B[N] by SUM.
          if ((alpha + 1.0) != 1.0)
            SUM = SUM * Gamma(alpha) * pow(x * 0.5, -alpha);
          TEMPA = ENMTEN;
          if (SUM > 1.0)
            TEMPA *= SUM;
          for (N = 0; N < nb; ++N) {
            if (fabs(B[N]) < TEMPA)
              B[N] = 0.0;
            B[N] = B[N] / SUM;
          }
          return NCALC;
        }
        M = 2 - M;
        if (M != 0) {
          EM -= 1.0;
          ALP2EM = (EM + EM) + alpha;
          ALPEM = (EM - 1.0) + alpha;
          if (ALPEM == 0)
            ALPEM = 1.0;
          SUM = (SUM + B[N - 1] * ALP2EM) * ALPEM / EM;
        }
      }
    }

    NEND = N - 2;

    if (NEND != 0) {
      // Calculate via difference equation and store B[N - 1], until N = 2
      for (L = 1; L <= NEND; ++L) {
        N -= 1;
        EN -= 2.0;
        B[N - 1] = (EN * B[N]) / x - B[N + 1];
        M = 2 - M;
        if (M != 0) {
          EM -= 1;
          ALP2EM = (EM + EM) + alpha;
          ALPEM = (EM - 1) + alpha;
          if (ALPEM == 0)
            ALPEM = 1.0;
          SUM = (SUM + B[N - 1] * ALP2EM) * ALPEM / EM;
        }
      }
    }

    // Calculate B[0]
    B[0] = 2.0 * (alpha + 1.0) * B[1] / x - B[2];

    EM -= 1; // statement 240
    ALP2EM = (EM + EM) + alpha;
    if (ALP2EM == 0)
      ALP2EM = 1;
    SUM += B[0] * ALP2EM;

    // Normalize. Divide all B[N] by SUM.
    if ((alpha + 1.0) != 1.0) // statement 250
      SUM = SUM * Gamma(alpha) * pow(x * 0.5, -alpha);
    TEMPA = ENMTEN;
    if (SUM > 1.0)
      TEMPA *= SUM;
    for (N = 0; N < nb; ++N) {
      if (fabs(B[N]) < TEMPA)
        B[N] = 0.0;
      B[N] = B[N] / SUM;
    }
  }

  return NCALC;
}

int bessel_Yn(int NB, double ALPHA, double X, double *BY) {
  int I, NA, NCALC;
  double ALFA, AYE, B, C, COSMU, D, DEN, DDIV, DIV, DMU, D1, D2,
    E, EN, ENU, EN1, EVEN, EX, F, G, GAMMA, H, ODD, P, PA, PA1,
    Q, QA, QA1, Q0, R, S, SINMU, TERM, TWOBYX, XNA, X2, YA, YA1;
  const double FIVPI = 1.5707963267948966192e1;  // 5 * pi
  const double ONBPI = 3.1830988618379067154e-1; // 1 / pi
  const double PI = 3.1415926535897932385e0;     // pi
  const double PIBY2 = 1.5707963267948966192e0;  // pi / 2
  const double PIM5 = 7.0796326794896619231e-1;  // 5 * pi - 15
  const double SQ2BPI = 7.9788456080286535588e-1; // sqrt(2 / pi)
  const double THRESH = 16.0;
  const double DEL = 1.0e-8;
  const double XMIN = 4.46e-308;
  const double XINF = 1.79e308;
  const double EPS = 1.11e-16;
  const double XLARGE = 1.0e8;

  // Coefficients for Chebyshev polynomial expansion of 1 / Gamma(1 - x)
  // where |x| <= 5
  const double CH[21] =
    {-0.67735241822398840964e-23,-0.61455180116049879894e-22,
     0.29017595056104745456e-20, 0.13639417919073099464e-18,
     0.23826220476859635824e-17,-0.90642907957550702534e-17,
     -0.14943667065169001769e-14,-0.33919078305362211264e-13,
     -0.17023776642512729175e-12, 0.91609750938768647911e-11,
     0.24230957900482704055e-09, 0.17451364971382984243e-08,
     -0.33126119768180852711e-07,-0.86592079961391259661e-06,
     -0.49717367041957398581e-05, 0.76309597585908126618e-04,
     0.12719271366545622927e-02, 0.17063050710955562222e-02,
     -0.76852840844786673690e-01,-0.28387654227602353814e+00,
     0.92187029365045265648e+00};

  EX = X;
  ENU = ALPHA;

  // Return if inputs are out of range
  if (NB <= 0 || X < XMIN || EX >= XLARGE || ENU < 0 || ENU >= 1) {
    return -1;
  }

  XNA = (int) (ENU + 0.5);
  NA = (int) (XNA);

  if (NA == 1)
    ENU = ENU - XNA;

  if (ENU == -0.5) {
    P = SQ2BPI / sqrt(EX);
    YA = P * sin(EX);
    YA1 = -P * cos(EX);
  } else if (EX < 3.0) {
    // Use Temme's scheme for small x
    B = EX * 0.5;
    D = -log(B);
    F = ENU * D;
    E = pow(B, -ENU);
    if (fabs(ENU) < DEL) {
      C = ONBPI;
    } else {
      C = ENU / sin(ENU * PI);
    }
    // Computation of sinh(f) / f
    if (fabs(F) < 1.0) {
      X2 = F * F;
      EN = 19.0;
      S = 1.0;
      for (I = 1; I <= 9; ++I) {
        S = S * X2 / EN / (EN - 1.0) + 1.0;
        EN = EN - 2.0;
      }
    } else {
      S = (E - 1.0 / E) * 0.5 / F;
    }

    // Computation of 1 / Gamma(1 - a) using Chebyshev polynomials
    X2 = ENU * ENU * 8.0;
    AYE = CH[0];
    EVEN = 0.0;
    ALFA = CH[1];
    ODD = 0.0;
    for (I = 3; I <= 19; I = I + 2) {
      EVEN = -(AYE + AYE + EVEN);
      AYE = -EVEN * X2 - AYE + CH[I - 1];
      ODD = -(ALFA + ALFA + ODD);
      ALFA = -ODD * X2 - ALFA + CH[I];
    }
    EVEN = (EVEN * 0.5 + AYE) * X2 - AYE + CH[20];
    ODD = (ODD + ALFA) * 2.0;
    GAMMA = ODD * ENU + EVEN;
    // End of computation of 1 / Gamma(1 - a)

    G = E * GAMMA;
    E = (E + 1.0 / E) * 0.5;
    F = 2.0 * C * (ODD * E + EVEN * S * D);
    E = ENU * ENU;
    P = G * C;
    Q = ONBPI / G;
    C = ENU * PIBY2;
    if (fabs(C) < DEL) {
      R = 1.0;
    } else {
      R = sin(C) / C;
    }
    R = PI * C * R * R;
    C = 1.0;
    D = -B * B;
    H = 0.0;
    YA = F + R * Q;
    YA1 = P;
    EN = 0.0;
    while (fabs(G / (1.0 + fabs(YA))) +
           fabs(H / (1.0 + fabs(YA1))) > EPS) {
      EN += 1.0;
      F = (F * EN + P + Q) / (EN * EN - E);
      C = C * D / EN;
      P = P / (EN - ENU);
      Q = Q / (EN + ENU);
      G = C * (F + R * Q);
      H = C * P - EN * G;
      YA += G;
      YA1 += H;
    }
    YA = -YA;
    YA1 = -YA1 / B;
  } else if (EX < THRESH) {
    // Use Temme's scheme for moderate x
    C = (0.5 - ENU) * (0.5 + ENU);
    B = EX + EX;
    E = (EX * ONBPI * cos(ENU * PI) / EPS);
    E = E * E;
    P = 1.0;
    Q = -EX;
    R = 1.0 + EX * EX;
    S = R;
    EN = 2.0;
    while (R * EN * EN < E) {
      EN1 = EN + 1.0;
      D = (EN - 1.0 + C / EN) / S;
      P = (EN + EN - P * D) / EN1;
      Q = (-B + Q * D) / EN1;
      S = P * P + Q * Q;
      R = R * S;
      EN = EN1;
    }
    F = P / S;
    P = F;
    G = -Q / S;
    Q = G;
    EN -= 1.0;
    
    // Initializing EN1 to 0.0 to remove some warning, without affecting the
    // code behavior from its original Fortran implementation. 
    // 
    // In the Fortran code, if the previous while loop is not entered, EN1 is
    // not initialized when first entering the following while loop. That means
    // R and S are computed from uninitilized values, and subsequently D, P, Q,
    // F, and G are junk. Notice that the line EN -= 1.0 inside the following
    // while loop corresponds to the goto line immediately before the while
    // loop. A closer look of the Fortran code suggests that the first time the
    // following while loop is entered, the only meaningful thing it does is to
    // initilize EN1 with EN. Then EN will be decremented to enter the while
    // loop again, from then on, EN1 keeps the value of EN and 1. 
    // 
    // Setting EN1 to 0.0 will remove the warning on R and S. Then later EN1 is
    // initialized with the desired value. 
    EN1 = 0.0; 
    while (EN > 0.0) {
      R = EN1 * (2.0 - P) - 2.0;
      S = B + EN1 * Q;
      D = (EN - 1.0 + C / EN) / (R * R + S * S);
      P = D * R;
      Q = D * S;
      E = F + 1.0;
      F = P * E - G * Q;
      G = Q * E + P * G;
      EN1 = EN;
      EN -= 1.0;
    }
    F = 1.0 + F;
    D = F * F + G * G;
    PA = F / D;
    QA = -G / D;
    D = ENU + 0.5 - P;
    Q = Q + EX;
    PA1 = (PA * Q - QA * D) / EX;
    QA1 = (QA * Q + PA * D) / EX;
    B = EX - PIBY2 * (ENU + 0.5);
    C = cos(B);
    S = sin(B);
    D = SQ2BPI / sqrt(EX);
    YA = D * (PA * S + QA * C);
    YA1 = D * (QA1 * S - PA1 * C);
  } else {
    // Use Campbell's asymptotic scheme
    NA = 0;
    D1 = (int) (EX / FIVPI);
    I = (int) (D1);
    DMU = ((EX - 15.0 * D1) - D1 * PIM5) - (ALPHA + 0.5) * PIBY2;
    if (I - 2 * (I / 2) == 0) {
      COSMU = cos(DMU);
      SINMU = sin(DMU);
    } else {
      COSMU = -cos(DMU);
      SINMU = -sin(DMU);
    }
    DDIV = 8.0 * EX;
    DMU = ALPHA;
    DEN = sqrt(EX);
    // The following is the original implementation, initializing YA 
    // when K is 1 is always executed, but somehow the compiler generates an
    // warning that YA might be used uninitialized. 
    // 
    // for (K = 1; K <= 2; ++K) {
    //   P = COSMU;
    //   COSMU = SINMU;
    //   SINMU = -P;
    //   D1 = (2.0 * DMU - 1.0) * (2.0 * DMU + 1.0);
    //   D2 = 0.0;
    //   DIV = DDIV;
    //   P = 0.0;
    //   Q = 0.0;
    //   Q0 = D1 / DIV;
    //   TERM = Q0;
    //   for (I = 2; I <= 20; ++I) {
    //     D2 = D2 + 8.0;
    //     D1 = D1 - D2;
    //     DIV = DIV + DDIV;
    //     TERM = -TERM * D1 / DIV;
    //     P = P + TERM;
    //     D2 = D2 + 8.0;
    //     D1 = D1 - D2;
    //     DIV = DIV + DDIV;
    //     TERM = TERM * D1 / DIV;
    //     Q = Q + TERM;
    //     if (fabs(TERM) <= EPS)
    //       break;
    //   }
    //   P = P + 1.0;
    //   Q = Q + Q0;
    //   if (K == 1) {
    //     YA = SQ2BPI * (P * COSMU - Q * SINMU) / DEN;
    //   } else {
    //     YA1 = SQ2BPI * (P * COSMU - Q * SINMU) / DEN;
    //   }
    //   DMU = DMU + 1.0;
    // }

    // Explicitly do k = 1
    P = COSMU;
    COSMU = SINMU;
    SINMU = -P;
    D1 = (2.0 * DMU - 1.0) * (2.0 * DMU + 1.0);
    D2 = 0.0;
    DIV = DDIV;
    P = 0.0;
    Q = 0.0;
    Q0 = D1 / DIV;
    TERM = Q0;
    for (I = 2; I <= 20; ++I) {
      D2 = D2 + 8.0;
      D1 = D1 - D2;
      DIV = DIV + DDIV;
      TERM = -TERM * D1 / DIV;
      P = P + TERM;
      D2 = D2 + 8.0;
      D1 = D1 - D2;
        DIV = DIV + DDIV;
        TERM = TERM * D1 / DIV;
        Q = Q + TERM;
        if (fabs(TERM) <= EPS)
          break;
    }
    P = P + 1.0;
    Q = Q + Q0;
    YA = SQ2BPI * (P * COSMU - Q * SINMU) / DEN;
    DMU = DMU + 1.0;
 
    // Explicitly do k = 2 
    P = COSMU;
    COSMU = SINMU;
    SINMU = -P;
    D1 = (2.0 * DMU - 1.0) * (2.0 * DMU + 1.0);
    D2 = 0.0;
    DIV = DDIV;
    P = 0.0;
    Q = 0.0;
    Q0 = D1 / DIV;
    TERM = Q0;
    for (I = 2; I <= 20; ++I) {
      D2 = D2 + 8.0;
      D1 = D1 - D2;
      DIV = DIV + DDIV;
      TERM = -TERM * D1 / DIV;
      P = P + TERM;
      D2 = D2 + 8.0;
      D1 = D1 - D2;
      DIV = DIV + DDIV;
      TERM = TERM * D1 / DIV;
      Q = Q + TERM;
      if (fabs(TERM) <= EPS)
        break;
    }
    P = P + 1.0;
    Q = Q + Q0;
    YA1 = SQ2BPI * (P * COSMU - Q * SINMU) / DEN;
    DMU = DMU + 1.0;
  }

  if (NA == 1) {
    H = 2.0 * (ENU + 1.0) / EX;
    if (H > 1.0) {
      if (fabs(YA1) > XINF / H) {
        H = 0.0;
        YA = 0.0;
      }
    }
    H = H * YA1 - YA;
    YA = YA1;
    YA1 = H;
  }

  // Now have first one or two Y's
  BY[0] = YA;
  BY[1] = YA1;
  if (YA1 == 0.0) {
    NCALC = 1;
  } else {
    AYE = 1.0 + ALPHA;
    TWOBYX = 2.0 / EX;
    NCALC = 2;
    for (I = 3; I <= NB; ++I) {
      if (TWOBYX < 1.0) {
        if (fabs(BY[I - 2]) * TWOBYX >= XINF / AYE) {
          break; // goto 450
        }
      } else {
        if (fabs(BY[I - 2]) >= XINF / AYE / TWOBYX) {
          break; // goto 450
        }
      }
      BY[I - 1] = TWOBYX * AYE * BY[I - 2] - BY[I - 3];
      AYE += 1.0;
      NCALC++;
    }
  }

  // statement 450
  for (I = NCALC + 1; I <= NB; ++I)
    BY[I - 1] = 0.0;

  return NCALC;
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

void bessel_jn_scaled(int nb, double x, double scale, double *B) {
  const double ensig = 1e-4;
  const double enmten = 1e-300;

  if (x <= ensig) {
    // Use 2-term Taylor expansion
    double scale_x = x / scale;
    double t1 = 1.0;
    double t2 = 0.5 * x * x;
    B[0] = t1 * (1 - t2 / 3.0);
    for (int i = 1; i <= nb; ++i) {
      t1 = t1 * scale_x / (2 * i + 1);
      if (t1 <= enmten) {
        t1 = 0.0;
      }
      B[i] = t1 * (1.0 - t2 / (2 * i + 3));
    }
  } else {
    double factor = sqrt(M_PI_2 / x);

    assert(bessel_Jn(nb + 1, 0.5, x, B) == (nb + 1));

    for (int i = 0; i <= nb; ++i) {
      B[i] *= factor;
      factor /= scale;
      if (fabs(B[i]) <= enmten) {
        factor = 0.0;
      }
    }
  }
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

void bessel_hn_scaled(int nb, double x, double scale, dcomplex_t *B) {
  const double ensig = 1e-4;
  const double enmten = 1e-300;
  double *Jn = new double[nb + 1];
  double *Yn = new double[nb + 1];

  if (x < ensig) {
    // User 2-term Taylor expansion
    double xscale1 =  x * scale;
    double xscale2 = scale / x;
    double t1 = 1.0;
    double t2 = 1.0;
    double factor = 0.5 * x * x;
    Jn[0] = t1 * (1.0 - factor / 3.0);
    Yn[0] = -(1.0 - factor) / x;
    B[0] = dcomplex_t{Jn[0], Yn[0]};
    for (int i = 1; i <= nb; ++i) {
      t1 = t1 * xscale1 / (2 * i + 1);
      t2 = t2 * xscale2 * (2 * i - 1);
      if (t1 <= enmten)
        t1 = 0.0;
      Jn[i] = t1 * (1.0 - factor / (2 * i + 3));
      Yn[i] = -t2 * (1.0 - factor / (1 - 2 * i));
      B[i] = dcomplex_t{Jn[i], Yn[i]};
    }
  } else {
    double factor = sqrt(M_PI_2 / x);

    assert(bessel_Jn(nb + 1, 0.5, x, Jn) == (nb + 1));
    assert(bessel_Yn(nb + 1, 0.5, x, Yn) == (nb + 1));

    for (int i = 0; i <= nb; ++i)
      B[i] = dcomplex_t{Jn[i], Yn[i]};

    for (int i = 0; i <= nb; ++i) {
      B[i] = B[i] * factor;
      factor = factor * scale;
    }
  }
  delete[] Jn;
  delete[] Yn;
}

} // namespace dashmm


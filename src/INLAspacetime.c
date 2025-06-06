#include "INLAspacetime.h"

double pclogrange(double lrange, double lam, int dim) {
  // return log of the PC-prior density for the log of the
  // practical range parameter in the Mat'ern model defined as
  //     r = (8 \nu) / \kappa
  // See Lindgren and Rue (2015) for this parametrization
  // See Fuglstad et. al. (2018) for this prior definition
  double dh = 0.5 * ((double)dim);
  return log(lam) -dh * lrange - lam * exp(-dh * lrange) + log(dh);
}

double pclogsigma(double lsigma, double lam) {
  // return log of the PC-prior density for the log of the
  // standard deviation parameter.
  // See Simpson et. al. (2017) for this prior definition
  return log(lam) + lsigma - lam * exp(lsigma);
}

void CSphere_gamma_alpha(double *lnGamma2, double *dalpha, double *cska) {
  // Compute the C(\gamma, \alpha) part depending on \gamma and \alpha
  // in eq. (23) of Lindgren et. al. (2024) <doi: 10.57645/20.8080.02.13>
  // \gamma and \alpha.
  //  1/(\gamma^(2\alpha))              if   \gamma^2 < exp(-0.7),
  //  1/((\alpha-1)\gamma^(2\alpha-2))  if   \gamma^2 > exp(0.7),
  //  sum_k (1+2k)/[\gamma^2 + k(k+1)]^\alpha, otherwise.
  // Note:  1/(4pi) is not included here.
  double out = 0.0;
  if((*lnGamma2) < (-1.4)) {
    out = -(*dalpha) * (*lnGamma2);
  } else {
    if((*lnGamma2) > 1.4) {
      out = -log((*dalpha)-1.0) -((*dalpha)-1.0)*(*lnGamma2);
    } else {
      double g2 = exp(*lnGamma2);			       // gamma^2
      int ialpha = ((int)(*dalpha));
      int isalphaInt = fabs(((double)ialpha) - (*dalpha)) < 0.0001;
      if (ialpha > 4L)
        isalphaInt = 0L;
      if (isalphaInt) {
        if (ialpha == 1L) {
          for (int k = 0; k < 50; k++) {
            out += (1 + 2*((double)k)) / (g2 + ((double)((k*(k+1)))));
          }
        }
        if (ialpha == 2L) {
          for (int k = 0; k < 50; k++) {
            out += (1 + 2*((double)k)) / pow2(g2 + ((double)((k*(k+1)))));
          }
        }
        if (ialpha == 3L) {
          for (int k = 0; k < 50; k++) {
            out += (1 + 2*((double)k)) / pow3(g2 + ((double)((k*(k+1)))));
          }
        }
        if (ialpha == 4L) {
          for (int k = 0; k < 50; k++) {
            out += (1 + 2*((double)k)) / pow4(g2 + ((double)((k*(k+1)))));
          }
        }
      } else {
        for (int k = 0; k < 50; k++) {
          out += (1 + 2*((double)k)) / pow(g2 + ((double)((k*(k+1)))), (*dalpha));
        }
      }
      out = log(out);
    }
  }
  cska[0] = out;
}

void ar2cov(int *n, int *k,
            double *a1, double *a2,
            double *r) {
  // function to compute the autocovariance for
  // the second order autoregression model - AR2
  int i, j, l0=2, l1=1, l2=0;
  double a, b, r1, r2;
  for(i=0; i<(*n); i++) {
    a = a1[i];
    b = a2[i];
    for(j=2; j<(*k); j++) {
      r2 = b * r[l2++];
      r1 = a * r[l1++];
      r[l0++] = r1 + r2;
    }
  }
}

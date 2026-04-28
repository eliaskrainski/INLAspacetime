#include "INLAspacetime.h"

double pclogrange(double lrange, double lam, int dim) {
  // return log of the PC-prior density for the log of the
  // practical range parameter in the Mat'ern model defined as
  //     r = (8 \nu) / \kappa
  // See Lindgren and Rue (2015) for this parametrization
  // See Fuglstad et. al. (2018) for this prior definition
  double dh = 0.5 * ((double)dim);
  return log(lam * dh) - dh * lrange - lam * exp(-dh * lrange);
}

double pclogsigma(double lsigma, double lam) {
  // return log of the PC-prior density for the log of the
  // standard deviation parameter.
  // See Simpson et. al. (2017) for this prior definition
  return log(lam) - lam * exp(lsigma) + lsigma;
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
      double den, g2 = exp(*lnGamma2);	// gamma^2
      int ialpha = ((int)(*dalpha));
      int isalphaInt = fabs(((double)ialpha) - (*dalpha)) < 0.0001;
      if (ialpha > 4L)
        isalphaInt = 0L;
      if (isalphaInt) {
        if (ialpha == 1L) {
          for (int k = 0; k < 50; k++) {
            den = (g2 + ((double)((k*(k+1)))));
            out += (1 + 2*((double)k)) / den;
          }
        }
        if (ialpha == 2L) {
          for (int k = 0; k < 50; k++) {
            den = pow2(g2 + ((double)((k*(k+1)))));
            out += (1 + 2*((double)k)) / den;
          }
        }
        if (ialpha == 3L) {
          for (int k = 0; k < 50; k++) {
            den = pow3(g2 + ((double)((k*(k+1)))));
            out += (1 + 2*((double)k)) / den;
          }
        }
        if (ialpha == 4L) {
          for (int k = 0; k < 50; k++) {
            den = pow4(g2 + ((double)((k*(k+1)))));
            out += (1 + 2*((double)k)) / den;
          }
        }
      } else {
        for (int k = 0; k < 50; k++) {
          den = pow(g2 + ((double)((k*(k+1)))), (*dalpha));
          out += (1 + 2*((double)k)) / den;
        }
      }
      out = log(out);
    }
  }
  cska[0] = out;
}

void ar2covk(int *n, int *k,
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
void cWMatern(int *N, double *S2, double *Scale,
              double *Nu, double *x, double *cc) {
  int n = (*N);
  double sigma2 = (*S2), nu=(*Nu), kappa = (*Scale);
  if(nu==0.5) {
    for(int i=0; i<n; i++) {
      cc[i] = sigma2 * exp(-kappa * x[i]);
    }
  } else {
    double hs, hsmall = 0.000000001;
    double scale, expon, cceps, bsmall;
    scale  = pow(2, 1.0 - nu) / gammafn(nu);
    double sig2scale = sigma2 * scale;
    cceps  = scale * pow(hsmall, nu) * bessel_k(hsmall, nu, 1);
    if(nu<1) {
      expon = 2.0*nu;
    } else {
      expon = 2.0;
    }
    bsmall = (1.0 - cceps) / pow(hsmall, expon);
    for(int i=0; i<n; i++) {
      hs = kappa * x[i];
      if(hs<hsmall) {
        if(hs==0) {
          cc[i] = sigma2;
        } else {
          cc[i] = sigma2 * (1.0 - bsmall * pow(hs, expon));
        }
      } else {
        cc[i] = sig2scale * pow(hs, nu) * bessel_k(hs, nu, 1);
      }
    }
  }
}

void c2ad(int *N, int *M, int *Mi, // inputs
          double *cb, double *cc,  // inputs
          double *d, double *aa) { // [in/ou]
  // cc contains all  C[nni, nni]
  // input cb = aa  = C[nni, i+1]
  // input d        = C[i+1, i+1]
  int i, j, l, m = *M, mi, k1=0, k2=0, one=1, info;
  double cw[m*m];
  char uplo = 'L'; // Fortran is column major order
  for(i=0; i<(*N); i++) {
    mi = Mi[i]; // Note: if i=0, mi=0
    if(mi>0) {
      for(j=0; j<mi; j++) {
        for(l=j; l<mi; l++) {
          cw[j*mi + l] = cc[k1++]; // copy to L
        }
      }
      dposv_(&uplo, &mi, &one, &cw[0], &mi,
             &aa[k2], &mi, &info, F_ONE);
      d[i] -= ddot_(&mi, &cb[k2], &one, &aa[k2], &one);
      k2 += mi;
    }
  }
}

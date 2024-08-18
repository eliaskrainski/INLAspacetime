#include "ar2cov.h"

// function to compute the autocovariance for
// the second order autoregression model - AR2
void ar2cov(int *n, int *k,
            double *a1, double *a2,
            double *r) {
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

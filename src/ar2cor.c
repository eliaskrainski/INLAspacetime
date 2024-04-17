#include "ar2cor.h"

void ar2cor(int *n, int *k,
            double *a1, double *a2,
            double *r) {
  int i,j,l=2*(*k)*(*n);
   for(i=0; i<(*n); i++) {
     for(j=0; j<(*k); j++) {
       r[l] = a1[i] * r[l-1] + a2[i]*r[l-2];
       l++;
     }
   }
}

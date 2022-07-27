#include <assert.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include <cgeneric.h>

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))
#define abs(x) sqrt((x)*(-x))

# ifdef __SUPPORT_SNAN__
#  define iszero(x) (fpclassify (x) == FP_ZERO)
# else
#  define iszero(x) (((__typeof (x)) (x)) == 0)
# endif

inla_cgeneric_func_tp inla_cgeneric_ar2ss_model;



#include <assert.h>
//#if !defined(__FreeBSD__)
//#include <malloc.h>
//#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include <cgeneric.h>

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))


# ifdef __SUPPORT_SNAN__
#  define iszero(x) (fpclassify (x) == FP_ZERO)
# else
#  define iszero(x) (((__typeof (x)) (x)) == 0)
# endif

#if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#else
typedef int fortran_charlen_t;
#endif
#define F_ONE ((fortran_charlen_t)1)

void dgemv_(char* trans, int* M, int* N, double* alpha, double* A,
           int* LDA, double* x, int* incx,
           double* beta, double* y, int* incy,
	   fortran_charlen_t);
double ddot_(int* N, double* DX, int* INCX, double* DY, int*INCY);

inla_cgeneric_func_tp inla_cgeneric_ar2ss_model;

inla_cgeneric_func_tp inla_cgeneric_st121_model;



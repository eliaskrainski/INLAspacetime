
/* cgeneric_defs.h
 * 
 * Copyright (C) 2022-2022 Elias Krainski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Elias Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#ifndef __CGENERIC_DEFS_H__
#define __CGENERIC_DEFS_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS
#include <stddef.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "cgeneric.h"
#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))
#define pow2(x) SQR(x)
#define pow3(x) (pow2(x)*(x))
#define pow4(x) (pow2(x)*pow2(x))
#ifdef __SUPPORT_SNAN__
#define iszero(x) (fpclassify(x) == FP_ZERO)
#else
#define iszero(x) (((__typeof(x))(x)) == 0)
#endif
#if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#else
typedef int fortran_charlen_t;
#endif
#define F_ONE ((fortran_charlen_t)1)

void dgemv_(char *trans, int *M, int *N, double *alpha, double *A, int *LDA, double *x,
	    int *incx, double *beta, double *y, int *incy, fortran_charlen_t);
double ddot_(int *N, double *DX, int *INCX, double *DY, int *INCY);
inla_cgeneric_func_tp inla_cgeneric_ar2ss_model;
inla_cgeneric_func_tp inla_cgeneric_st121_model;
inla_cgeneric_func_tp inla_cgeneric_sstspde;

__END_DECLS
#endif

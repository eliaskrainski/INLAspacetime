
/* cgeneric_barrier.c
 *
 * Copyright (C) 2022-2023 Elias Krainski
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

#include "INLAspacetime.h"

// This function uses the padded matrices with zeroes

double *inla_cgeneric_barrier(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{

	double *ret = NULL;
	int N, M, ith, nth, ifix[2];
	double sigma2, range2, params[4];

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

/*
	assert(!strcasecmp(data->ints[1]->name, "debug"));     // this will always be the case
	int debug = data->ints[1]->ints[0];		       // this will always be the case
	assert(debug >= 0);				       // just to 'find an use for "debug"
*/

	assert(!strcasecmp(data->ints[2]->name, "ii"));
	inla_cgeneric_vec_tp *ii = data->ints[2];
	M = ii->len;

	assert(!strcasecmp(data->ints[3]->name, "jj"));
	inla_cgeneric_vec_tp *jj = data->ints[3];
	assert(M == jj->len);

	// prior parameters
	assert(!strcasecmp(data->doubles[0]->name, "prange"));
	inla_cgeneric_vec_tp *prange = data->doubles[0];
	assert(prange->len == 2);

	assert(!strcasecmp(data->doubles[1]->name, "psigma"));
	inla_cgeneric_vec_tp *psigma = data->doubles[1];
	assert(psigma->len == 2);

	nth = 0;
	if (iszero(prange->doubles[1])) {
		ifix[0] = 1;
	} else {
		ifix[0] = 0;
		nth++;
	}

	if (iszero(psigma->doubles[1])) {
		ifix[1] = 1;
	} else {
		ifix[1] = 0;
		nth++;
	}
	assert(nth < 3);

	if (theta) {
		// interpretable parameters
		// theta = log(range, sigma)
		ith = 0;
		if (ifix[0] == 1) {
		  range2 = SQR(prange->doubles[0]);
		} else {
			range2 = exp(theta[ith++] * 2.0);
		}

		if (ifix[1] == 1) {
		  sigma2 = SQR(psigma->doubles[0]);
		} else {
			sigma2 = exp(theta[ith++] * 2.0);
		}

		double pi2s2 = 0.6366197724 / sigma2;  // 2 / ( pi * sigma^2)

		params[0] = pi2s2 / (range2);
		params[1] = pi2s2 / (8);
		params[2] = pi2s2 / (8);
		params[3] = range2 * pi2s2 / (64);

		assert(nth == ith);

//		if (debug) {
//			fprintf(stderr, "theta = ");
//			for (int i = 0; i < nth; i++)
//				fprintf(stderr, "%f ", theta[i]);
//			fprintf(stderr, "range = %f, sigma = %f\n", range, sigma);
//		}
	} else {
		for (int i = 0; i < 4; i++) {
			params[i] = NAN;
		}
	}


	switch (cmd) {

	case INLA_CGENERIC_GRAPH:
	{
		int k = 2;
		ret = Calloc(k + 2 * M, double);
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (int i = 0; i < M; i++) {
			ret[k++] = ii->ints[i];
		}
		for (int i = 0; i < M; i++) {
			ret[k++] = jj->ints[i];
		}
	}
		break;

	case INLA_CGENERIC_Q:
	{
		int offset = 2;
		ret = Calloc(offset + M, double);

		assert(!strcasecmp(data->mats[0]->name, "xx"));
		inla_cgeneric_mat_tp *xx = data->mats[0];
		assert(xx->nrow == 4);
		assert(xx->ncol == M);

		int one = 1, nm = 4;
		double zerof = 0.0, onef = 1.0;
		char trans = 'N';

		dgemv_(&trans, &M, &nm, &onef, &xx->x[0], &M, &params[0], &one, &zerof, &ret[offset], &one, F_ONE);

		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;				       /* REQUIRED */
	}
		break;

	case INLA_CGENERIC_MU:
	{
		// return (N, mu). if N==0 then mu is not needed as its taken to be mu[]==0
		ret = Calloc(1, double);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(P, initials)
		// where P is the number of hyperparameters
		ret = Calloc(nth + 1, double);
		ith = 0;
		ret[ith++] = (double) nth;
		if (ifix[0] == 0) {
			ret[ith++] = 1.0;
		}
		if (ifix[1] == 0) {
			ret[ith++] = 1.0;
		}
		assert(ith == (nth + 1));
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
		// PC-priors
		ret[0] = 0.0;
		ith = 0;
		double lam; //, daux = 2 * 0.5;
		if (ifix[0] == 0) {
			lam = -log(prange->doubles[1]) / prange->doubles[0];
//			ret[0] += log(lam1 * daux) - daux * theta[ith] - lam1 * exp(-daux * theta[ith]) + log(daux);
			ret[0] += pclogrange(theta[ith], lam, 2L);
			ith++;
		}
		if (ifix[1] == 0) {
			lam = -log(psigma->doubles[1]) / psigma->doubles[0];
		  ret[0] += pclogsigma(theta[ith], lam);
		  ith++;
		}
		assert(ith == nth);
	}
		break;

	case INLA_CGENERIC_VOID:
	case INLA_CGENERIC_LOG_NORM_CONST:
	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}

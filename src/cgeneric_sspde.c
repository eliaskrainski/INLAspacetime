
/* cgeneric_sspde.c
 *
 * Copyright (C) 2025 Elias Krainski
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
double *inla_cgeneric_sspde(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{

	double *ret = NULL;
	int N, M, ith, nth, ifix[2];
	double lkappa2, ltau2;

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	assert(!strcasecmp(data->ints[1]->name, "debug"));     // this will always be the case
	int debug = data->ints[1]->ints[0];		       // this will always be the case
	assert(debug >= 0);				       // just to 'find an use for "debug" ...'
	if (debug>0) debug = 1;

	assert(!strcasecmp(data->ints[2]->name, "Rmanifold"));
	int Rmanifold = data->ints[2]->ints[0];
	assert(Rmanifold >= 0);

	assert(!strcasecmp(data->ints[3]->name, "dimension"));
	int dimension = data->ints[3]->ints[0];
	assert(dimension > 0);

	assert(!strcasecmp(data->ints[4]->name, "alpha"));
	int alpha = data->ints[4]->ints[0];
	double dalpha = (double)alpha;
	double nu_s = dalpha -0.5*((double)dimension);

	assert(!strcasecmp(data->ints[5]->name, "nm"));
	int nm = data->ints[5]->ints[0];
	assert(nm > 0);
	double params[nm];

	assert(!strcasecmp(data->ints[6]->name, "ii"));
	inla_cgeneric_vec_tp *ii = data->ints[6];
	M = ii->len;

	assert(!strcasecmp(data->ints[7]->name, "jj"));
	inla_cgeneric_vec_tp *jj = data->ints[7];
	assert(M == jj->len);

	assert(!strcasecmp(data->doubles[0]->name, "cc"));
	inla_cgeneric_vec_tp *cc = data->doubles[0];
	assert(cc->len == 2);

	// prior parameters for range
	assert(!strcasecmp(data->doubles[1]->name, "prange"));
	inla_cgeneric_vec_tp *prange = data->doubles[1];
	assert(prange->len == 2);

	// prior parameters for sigma
	assert(!strcasecmp(data->doubles[2]->name, "psigma"));
	inla_cgeneric_vec_tp *psigma = data->doubles[2];
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

//	FILE *fp = fopen("cg_stspde.log", "w");

	if (theta) {
	  // interpretable parameters: (range, sigma)
		// internal parameters: (log(kappa), log(tau))
		// as in Lindgren and Rue (2015), this is
		// the JSS paper https://doi.org/10.18637/jss.v063.i19
		ith = 0;
		if (ifix[0] == 1) {
		  lkappa2 = cc->doubles[0] - 2.0*log(prange->doubles[0]);
		} else {
		  lkappa2 = cc->doubles[0] - 2.0*theta[ith++];
		}

		// map to gamma_e depends on the manifold
		if (Rmanifold) { // if R_d (not sphere)
		  ltau2 = -nu_s * lkappa2;
		} else {// if sphere, add parameter dependent part of
		  // C_sphere(gamma_s=\kappa,alpha) in Eq. (23) of
		  // Lindgren et. al. (2024) <doi: 10.57645/20.8080.02.13>
		  CSphere_gamma_alpha(&lkappa2, &dalpha, &ltau2);
		}
		if (ifix[1] == 1) {
		  ltau2 += cc->doubles[1] - 2.0*log(psigma->doubles[0]);
		} else {
		  ltau2 += cc->doubles[1] - 2.0*theta[ith++];
		}
		assert(nth == ith);

	  params[2] = exp(ltau2);
	  params[1] = exp(ltau2 + lkappa2) * 2.0;
	  params[0] = exp(ltau2 + lkappa2 * 2.0);

	} else {
		for (int i = 0; i < nm; i++) {
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
		ret[0] = -1;				       /* REQUIRED */
    ret[1] = M;				       /* REQUIRED */

		assert(!strcasecmp(data->mats[0]->name, "xx"));
		inla_cgeneric_mat_tp *xx = data->mats[0];
		assert(xx->nrow == 3);
		assert(xx->ncol == M);

		int one = 1;
		double zerof = 0.0, onef = 1.0;
		char trans = 'N';

		dgemv_(&trans, &M, &nm, &onef, &xx->x[0], &M, params, &one, &zerof, &ret[offset], &one, F_ONE);

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
		double daux = 0.5 * ((double) dimension), lam;
		if (ifix[0] == 0) {
			lam = -log(prange->doubles[1]) * pow(prange->doubles[0], daux);
			ret[0] += pclogrange(theta[ith], lam, dimension);
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


/* cgeneric_ar2ss.c
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

#include <stdio.h>
#include "cgeneric_defs.h"

double *inla_cgeneric_ar2ss_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{

	double *ret = NULL, a1, a2, sleng, lprec, prec, pcorrect;

	assert(!strcasecmp(data->ints[0]->name, "n"));
	int N = data->ints[0]->ints[0];
	assert(N > 0);

	assert(!strcasecmp(data->doubles[0]->name, "psigma"));
	inla_cgeneric_vec_tp *psigma = data->doubles[0];
	assert(psigma->len == 2);

	assert(!strcasecmp(data->doubles[1]->name, "pleng"));
	inla_cgeneric_vec_tp *pleng = data->doubles[1];
	assert(pleng->len == 2);

	assert(!strcasecmp(data->doubles[2]->name, "pcor"));
	inla_cgeneric_vec_tp *pcor = data->doubles[2];
	assert(pcor->len == 2);

	assert(!strcasecmp(data->ints[1]->name, "toprint"));
	int toprint = data->ints[1]->ints[0];
	assert(toprint>=0);
	if(toprint>0) toprint = 1;

	int npar = 0;
	if (theta) {
		if (iszero(psigma->doubles[1])) {
			prec = 1 / SQR(psigma->doubles[0]);
			lprec = log(prec);
		} else {
			lprec = -2 * theta[npar];
			prec = exp(lprec);
			npar++;
		}
		if (iszero(pleng->doubles[1])) {
			sleng = pleng->doubles[0];
		} else {
			sleng = exp(theta[npar]);
			npar++;
		}
		if (iszero(pcor->doubles[1])) {
			a2 = pcor->doubles[0];
		} else {
			a2 = (2.0 / (1.0 + exp(-theta[npar]))) - 1.0;
			npar++;
		}
		a1 = -2.0 * sqrt(a2) * cos(2.0 * M_PI / sleng);
		pcorrect = exp(log(1 + a2) - log(1.0 - a2) - log(1.0 - a1 + a2) - log(1.0 + a1 + a2));
	} else {
		sleng = NAN;
		prec = lprec = NAN;
		a2 = NAN;
		pcorrect = a1 = NAN;
	}

	switch (cmd) {

	case INLA_CGENERIC_GRAPH:
	{
		int M = N + N - 1, k = 0;
		if (N > 2) {
			M += N - 2;
		}
		ret = Calloc(2 + 2 * M, double);

		ret[k++] = N;				       /* dimension */
		ret[k++] = M;				       /* number of (i <= j) */
		// printf("graph : N = %g, M = %g, k = %d\n", ret[0], ret[1], k);

		if (N == 1) {
			ret[k++] = 0;
			ret[k++] = 0;
		}
		if (N == 2) {
			ret[k++] = 0;
			ret[k++] = 0;
			ret[k++] = 1;
			ret[k++] = 0;
			ret[k++] = 1;
			ret[k++] = 1;
		}
		// printf("G k = %d\n", k);
		if (N > 2) {
			//if (toprint)
			//	fprintf(stderr, "graph\n");
			for (int i = 0; i < (N - 2); i++) {
				ret[k++] = i;		       /* i */
				ret[k++] = i;		       /* i */
				ret[k++] = i;		       /* i */
			}
			// printf("G k = %d\n", k);
			ret[k++] = N - 2;
			ret[k++] = N - 2;
			ret[k++] = N - 1;
			for (int i = 0; i < (N - 2); i++) {
				ret[k++] = i;		       /* j */
				ret[k++] = i + 1;	       /* j */
				ret[k++] = i + 2;	       /* j */
//				if (toprint)
//					fprintf(stderr, "%d %d %d %d\n", i, i, i + 1, i + 2);
			}
			ret[k++] = N - 2;
			ret[k++] = N - 1;
			ret[k++] = N - 1;
			// printf("G k = %d\n", k);
		}
		assert((2 * M + 2) == k);
	}
		break;

	case INLA_CGENERIC_Q:
	{
		double q0, q1, q2, param = pcorrect * prec;
		int M = N + N - 1;
		if (N > 2)
			M += N - 2;
		// printf("Q M: %d\n", M);
		ret = Calloc(2 + M, double);

		int k = 0;
		ret[k++] = -1;
		ret[k++] = M;
		// printf("Q : ret = %g, M = %g, k = %d\n", ret[0], ret[1], k);

//		if (toprint)
//			fprintf(stderr, "Q\n");

		if (N == 1) {
			ret[k++] = param;
		}
		if (N == 2) {
			ret[k++] = param;
			ret[k++] = 0.5 * a1 * param;	       // convention: 0.5 * rho * (1-rho^2)/s2innovation
			ret[k++] = param;
		}
		if (N > 2) {
			ret[k++] = param;
			ret[k++] = a1 * param;
			q2 = a2 * param;
			ret[k++] = q2;
			if (N > 3) {
				ret[k++] = (1.0 + a1 * a1) * param;
				q1 = a1 * (1.0 + a2) * param;
				ret[k++] = q1;
				ret[k++] = q2;
				if (N > 4) {
					q0 = (1.0 + a1 * a1 + a2 * a2) * param;
					for (int i = 3; i < (N - 1); i++) {
						ret[k++] = q0;
						ret[k++] = q1;
						ret[k++] = q2;
//						if (toprint)
//							fprintf(stderr, "%d %f %f %f\n", i, q0, q1, q2);
					}
				}
			}
			ret[k++] = (1.0 + a1 * a1) * param;
			ret[k++] = a1 * param;
			ret[k++] = param;
		}
	}
		;
		break;

	case INLA_CGENERIC_MU:
	{
		ret = Calloc(1, double);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		ret = Calloc(npar + 1, double);
		npar = 0;
		if (!iszero(psigma->doubles[1])) {
			ret[++npar] = 0.0;		       // log(psigma->doubles[0] - 1);
		}
		if (!iszero(pleng->doubles[1])) {
			ret[++npar] = 1.0;		       // log(pleng->doubles[0] + 1);
		}
		if (!iszero(pcor->doubles[1])) {
			ret[++npar] = 2.94445;		       // log((0.5+0.5*pcor->doubles[0])/(1-(0.5+0.5*pcor->doubles[0])));
		}
		ret[0] = npar;
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
		ret[0] = 0.0;
		double lam = 0.0;
		npar = 0;
		if (!iszero(psigma->doubles[1])) {
			lam = -log(psigma->doubles[1]) / psigma->doubles[0];
			ret[0] += log(lam) + theta[npar] - lam * exp(theta[npar]);
			npar++;
		}
		if (!iszero(pleng->doubles[1])) {
			lam = -log(pleng->doubles[1]) / pleng->doubles[0];
			ret[0] += log(lam * 0.5) - 0.5 * theta[npar] - lam * exp(-0.5 * theta[npar]);
			npar++;
		}
		if (!iszero(pcor->doubles[1])) {
			// ret[0] += priorfunc_pc_cor1(&theta[npar], &pcor->doubles[0]);
			ret[0] += -0.5 * SQR(theta[npar]);
			npar++;
		}
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

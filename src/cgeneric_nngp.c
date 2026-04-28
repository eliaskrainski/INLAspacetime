
/* cgeneric_nngp.c
 *
 * Copyright (C) 2026-2027 Elias Krainski
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

double *inla_cgeneric_nngp(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{

	double *ret = NULL;
	int i, j, N, M, ith, nth, ifix[2];

	assert(data->n_ints > 1);

	// the size of the model
	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	assert(!strcasecmp(data->ints[1]->name, "debug"));     // this will always be the case
	int debug = data->ints[1]->ints[0];		       // this will always be the case
	assert(debug >= 0);			// just to 'find an use for "debug" ...'
	if (debug>0) debug = 1; // just to 'find an use for "debug" ...'

	// correlation function
	assert(!strcasecmp(data->ints[2]->name, "cfn"));
	int cfn = data->ints[2]->ints[0];

	assert(!strcasecmp(data->ints[3]->name, "ii"));
	inla_cgeneric_vec_tp *ii = data->ints[3];
	M = ii->len;

	assert(!strcasecmp(data->ints[4]->name, "jj"));
	inla_cgeneric_vec_tp *jj = data->ints[4];
	assert(M == jj->len);

	// prior parameters for range
	assert(!strcasecmp(data->doubles[0]->name, "prange"));
	inla_cgeneric_vec_tp *prange = data->doubles[0];
	assert(prange->len == 2);

	// prior parameters for sigma
	assert(!strcasecmp(data->doubles[1]->name, "psigma"));
	inla_cgeneric_vec_tp *psigma = data->doubles[1];
	assert(psigma->len == 2);

	// smoothness nu (fixed)
	assert(!strcasecmp(data->doubles[2]->name, "nu"));
	double nu = data->doubles[2]->doubles[0];


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

	// phi = 1/kappa;
	// range = sqrt(8 * nu) / kappa
	// phi = sqrt(8 * nu) * range
	double scale, sigma2;
	if (theta) {
		ith = 0;
	  if(cfn==1) {
	    if (ifix[0] == 1) {
	      scale = sqrt(8.0 * nu) / prange->doubles[0];
	    } else {
	      scale = sqrt(8.0 * nu) / exp(theta[ith++]);
	    }
	  } else {
	    if (ifix[0] == 1) {
	      // here there is close expression for range at correl = 0.05
	      // but here we use 0.135 because this is the case for Matern
	      scale = pow(2.0, 1/nu) / prange->doubles[0];
	    } else {
	      scale = pow(2.0, 1/nu) / exp(theta[ith++]);
	    }
	  }
		if (ifix[1] == 1) {
		  sigma2 = SQR(psigma->doubles[0]);
		} else {
		  sigma2 = exp(2 * theta[ith++]);
		}
		assert(nth == ith);

	} else {
	  scale = NAN;
	  sigma2 = NAN;
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
		ret[0] = -1;	// REQUIRED
    ret[1] = M;		// REQUIRED
    double daux;

    assert(!strcasecmp(data->ints[5]->name, "Mmax"));
    int Mmax = data->ints[5]->ints[0];
    assert(!strcasecmp(data->ints[6]->name, "Mi"));
    inla_cgeneric_vec_tp *Mi = data->ints[6];
    assert(!strcasecmp(data->ints[7]->name, "Aj"));
    inla_cgeneric_vec_tp *Aj = data->ints[7];
    assert(!strcasecmp(data->ints[8]->name, "nll"));
    inla_cgeneric_vec_tp *nll = data->ints[8];
    assert(!strcasecmp(data->ints[9]->name, "iL2"));
    inla_cgeneric_vec_tp *iL1 = data->ints[9];
    assert(!strcasecmp(data->ints[10]->name, "iL1"));
    inla_cgeneric_vec_tp *iL2 = data->ints[10];

    assert(!strcasecmp(data->doubles[3]->name, "cbdists"));
//    double *cbd = &data->doubles[3]->doubles[0];
    int ncb = data->doubles[3]->len;
    double cb[ncb], aa[ncb];
    assert(!strcasecmp(data->doubles[4]->name, "ccdists"));
  //  double *ccd = &data->doubles[4]->doubles[0];
    int ncc = data->doubles[4]->len;
    double d[N], cc[ncc];
    if(cfn==1) {
      cWMatern(&ncb, &sigma2, &scale, &nu,
               &data->doubles[3]->doubles[0], &cb[0]);
      cWMatern(&ncc, &sigma2, &scale, &nu,
               &data->doubles[4]->doubles[0], &cc[0]);
    } else {
      if(nu==1) {
        for(i=0; i<ncb; i++) {
          cb[i] = sigma2 * exp(-data->doubles[3]->doubles[i]*scale);
        }
        for(i=0; i<ncc; i++) {
          cc[i] = sigma2 * exp(-data->doubles[4]->doubles[i]*scale);
        }
      } else {
        for(i=0; i<ncb; i++) {
          daux = pow(data->doubles[3]->doubles[i]*scale, nu);
          cb[i] = sigma2 * exp(-daux);
        }
        for(i=0; i<ncc; i++) {
          daux = pow(data->doubles[4]->doubles[i]*scale, nu);
          cc[i] = sigma2 * exp(-daux);
        }
      }
    }
    for(i=0; i<ncb; i++) {
      aa[i] = cb[i];
    }
    for(i=0; i<N; i++) {
      d[i] = sigma2;
    }
    // compute A and D
    c2ad(&N, &Mmax, &Mi->ints[0],
         &cb[0], &cc[0], &d[0], &aa[0]);
    // compute L
    double ll[N+ncb], sdi;
    int k1=0, k2=0; j=0;
    for(i=0; i<N; i++) {
      sdi = sqrt(d[i]); //(Aj->ints[k1]]);
      if(Mi->ints[i]>0) {
        for(j=0; j<Mi->ints[i]; j++) {
          ll[k2++] = -aa[k1++] / sdi;
        }
      }
      ll[k2++] = 1/sdi; // diag of ll, the last at each line
    }
/*    FILE *fp = fopen("ll.log", "w");
    k1=0; k2=0;
    for(i=0; i<N; i++) {
      if(Mi->ints[i]>0) {
        for(j=0; j<Mi->ints[i]; j++) {
          fprintf(fp, "%d %d %f\n", i, Aj->ints[k1++], ll[k2++]);
        }
      }
      fprintf(fp, "%d %d %f\n", i, i, ll[k2++]);
    }
    fclose(fp);
*/
    // compute LL'
    int k = 0;
    for(i=0; i<M; i++) { // each element of Q, Q[i,j]
      daux = 0;
      for(j=0; j<nll->ints[i]; j++) {
        daux += ll[iL1->ints[k]] * ll[iL2->ints[k]];
        k++;
      }
      ret[offset+i] = daux;
    }
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
			ret[ith++] = 10.0;
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
		int dimension = 2;
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

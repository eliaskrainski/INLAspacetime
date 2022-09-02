#include "cgeneric_defs.h"

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_st121_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  // define model121 constants
  //  double at = 1, as=2, ae=1, aa = 2, nus = 1;
  // c1 = gamma(a_t - 1/2) / [ gamma(a_t)(4 pi)^(1/2) ]
  // c2 = gamma(a - 1) / [ gamma(a)(4pi) ]
  double lc12 = -3.224171; //log(0.03978874); //4 * 4.141592654;
  double lg[3], params[9];
  int N, M, nc, i;

  // the size of the model
  assert(data->n_ints > 1);
  assert(!strcasecmp(data->ints[0]->name, "n"));       // this will always be the case
  N = data->ints[0]->ints[0];			       // this will always be the case
  assert(N > 0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));    // this will always be the case
  int debug = data->ints[1]->ints[0];	        // this will always be the case

  if(debug==1) {
    debug = 1; // just to 'find an use for "debug" ...'
  }

  assert(!strcasecmp(data->ints[2]->name, "ii"));
  inla_cgeneric_vec_tp *ii = data->ints[2];
  M = ii->len;

  assert(!strcasecmp(data->ints[3]->name, "jj"));
  inla_cgeneric_vec_tp *jj = data->ints[3];
  assert(M == jj->len);

  assert(!strcasecmp(data->mats[0]->name, "xx"));
  inla_cgeneric_mat_tp *xx = data->mats[0];
  assert(M == xx->ncol);
  nc = xx->nrow;

  // prior parameters
  assert(!strcasecmp(data->doubles[0]->name, "prs"));
  inla_cgeneric_vec_tp *prs = data->doubles[0];
  assert(prs->len == 2);

  assert(!strcasecmp(data->doubles[1]->name, "prt"));
  inla_cgeneric_vec_tp *prt = data->doubles[1];
  assert(prt->len == 2);

  assert(!strcasecmp(data->doubles[2]->name, "psigma"));
  inla_cgeneric_vec_tp *psigma = data->doubles[2];
  assert(psigma->len == 2);

  int ith, ifix[3];
  if(iszero(prs->doubles[1])){
    ifix[0] = 1;
  } else {
    ifix[0] = 0;
  }
  if(iszero(prt->doubles[1])) {
    ifix[1] = 1;
  } else {
    ifix[1] = 0;
  }
  if(iszero(psigma->doubles[1])) {
    ifix[2] = 1;
  } else {
    ifix[2] = 0;
  }
  int nth = 3;
  for(i=0; i<3; i++)
    nth -= ifix[i];

  if (theta) {
    // interpretable parameters
    //  theta = log(range_s, range_t, sigma)
    //   lg = log(g_s, g_t, g_e)
    // g_s = sqrt(8 * v_s) / r_s ;
    // g_t = r_t * g_s^a_s / sqrt(8(a_t-1/2)) ;
    // g_e = sqrt(c12 / (sigma * g_t * g_s^{2a-d})) ;
    ith = 0;
    if(ifix[0]){
      lg[0] = 1.039721 - log(prs->doubles[0]);
    } else {
      lg[0] = 1.039721 - theta[ith];
      ith++;
    }
    if(iszero(prt->doubles[1])) {
      lg[1] = log(prt->doubles[0]) + lg[0]*2 - 0.6931472;
    } else {
      lg[1] = theta[ith] + lg[0]*2 - 0.6931472;
      ith++;
    }
    if(iszero(psigma->doubles[1])) {
      lg[2] = 0.5*( lc12 - log(psigma->doubles[0]) - lg[0] - lg[1]*2 ) ;
    } else {
      lg[2] = 0.5*( lc12 - theta[ith]*2 - lg[0] - lg[1]*2 ) ;
      ith++;
    }
    assert(nth == ith);

    params[0] = exp(lg[2]*2 + lg[0]*6);             //  g_s^6         * g_e^2
    params[1] = exp(lg[2]*2 + lg[0]*4 + 3);         // 3g_s^4         * g_e^2
    params[2] = exp(lg[2]*2 + lg[0]*2 + 3);         // 3g_s^2         * g_e^2
    params[3] = exp(lg[2]*2);                       //      1         * g_e^2
    params[4] = exp(lg[2]*2 + lg[1] + lg[0]*4);     //  g_s^4 * g_t   * g_e^2
    params[5] = exp(lg[2]*2 + lg[1] + lg[0]*2 + 2); // 2g_s^2 * g_t   * g_e^2
    params[6] = exp(lg[2]*2 + lg[1]);               //          g_t   * g_e^2
    params[7] = exp(lg[2]*2 + lg[1]*2 + lg[0]*2);   //  g_s^2 * g_t^2 * g_e^2
    params[8] = exp(lg[2]*2 + lg[1]*2);             //          g_t^2 * g_e^2
  } else {
    params[3] = params[2] = params[1] = params[0] = NAN;
    params[6] = params[5] = params[4]  = NAN;
    params[8] = params[7] = NAN;
  }


  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }

  case INLA_CGENERIC_GRAPH:
    {
      int offset = 2;
      ret = Calloc(2 + 2 * M, double);
      ret[0] = N;        	       /* dimension */
      ret[1] = M;		   /* number of (i <= j) */
      for (i = 0; i < M; i++) {
	ret[offset + i] = ii->ints[i];
	ret[offset + M + i] = jj->ints[i];
      }
      break;
    }

  case INLA_CGENERIC_Q:
    {
      int offset = 2;
      ret = Calloc(2 + M, double);
      //memset(ret + offset, 0, M * sizeof(double));
      int one=1;
      double zerof = 0.0, onef = 1.0;
      char trans  = 'N'; // (trans=="N") | (trans=='n')) {
      dgemv_(&trans, &M, &nc, &onef, &xx->x[0], &M, params, &one, &zerof, &ret[offset], &one, F_ONE);
      ret[0] = -1;		/* REQUIRED */
      ret[1] = M;		/* REQUIRED */
      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (N, mu)
      // if N==0 then mu is not needed as its taken to be mu[]==0
      ret = Calloc(1, double);
      ret[0] = 0;
      break;
    }

  case INLA_CGENERIC_INITIAL:
    {
      // return c(M, initials)
      // where M is the number of hyperparameters
      ret = Calloc(nth+1, double);
      ret[0] = (double)nth;
      ith = 0;
      if(ifix[0]==0){
	ith++;
	ret[ith] = log(prs->doubles[0]) + 1;
      }
      if(ifix[1]==0) {
	ith++;
	ret[ith] = log(prt->doubles[0]) + 1;
      }
      if(ifix[2]==0) {
	ith++;
	ret[ith] = log(psigma->doubles[0]) -1.0;
      }
      assert(ith==nth);
      break;
    }

  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      ret = Calloc(1, double);
      // PC-priors
      ret[0] = 0.0;
      ith = 0;
      if(ifix[0]==0) {
	double lam2 = -log(prs->doubles[1])/prs->doubles[0];
	ret[0] += log(lam2) - theta[ith] - lam2*exp(-1*theta[ith]) ;
	ith++;
      }
      if(ifix[1]==0) {
	double lam3 = -log(prt->doubles[1])/prt->doubles[0];
	ret[0] += log(lam3) - 0.5*theta[ith] - lam3*exp(-0.5*theta[ith]) - log(0.5) ;
	ith++;
      }
      if(ifix[2]==0) {
	double lam1 = -log(psigma->doubles[1])/psigma->doubles[0];
	ret[0] += log(lam1) + theta[ith] - lam1*exp(theta[ith]) ;
	ith++;
      }
      assert(ith==nth);
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}


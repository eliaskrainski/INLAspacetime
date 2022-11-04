#include "cgeneric_defs.h"

// This function uses the padded matrices with zeroes
double *inla_cgeneric_sstspde(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  int N, M, i, k, ith, nth, ifix[3];
  double lg[3], a1, a2;
  
  // the size of the model
  assert(data->n_ints > 1);
  assert(!strcasecmp(data->ints[0]->name, "n"));       // this will always be the case
  N = data->ints[0]->ints[0];			       // this will always be the case
  assert(N > 0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));    // this will always be the case
  int debug = data->ints[1]->ints[0];	        // this will always be the case
  assert(debug>=0); // just to 'find an use for "debug" ...'
    
  assert(!strcasecmp(data->ints[2]->name, "ii"));
  inla_cgeneric_vec_tp *ii = data->ints[2];
  M = ii->len;

  assert(!strcasecmp(data->ints[3]->name, "jj"));
  inla_cgeneric_vec_tp *jj = data->ints[3];
  assert(M == jj->len);

  assert(!strcasecmp(data->ints[4]->name, "aaa"));
  inla_cgeneric_vec_tp *aaa = data->ints[4];
  assert(aaa->len == 3);
  double alphas = (float)aaa->ints[1];
  double alpha = (float)aaa->ints[2] + alphas*((float)aaa->ints[0]-0.5);
  int ialpha = ((float)((int)alpha)) == alpha;
  double aaux = 1-alpha; // d=2;

  assert(!strcasecmp(data->ints[5]->name, "manifold"));
  int manifold = data->ints[5]->ints[0];
  assert(manifold>=0);

  assert(!strcasecmp(data->ints[6]->name, "nm"));
  int nm = data->ints[6]->ints[0];
  assert(nm>0);
  double params[nm];

  assert(!strcasecmp(data->doubles[0]->name, "cc"));
  inla_cgeneric_vec_tp *cc = data->doubles[0];
  assert(cc->len == 3);
  double c3 = cc->doubles[2];

  assert(!strcasecmp(data->doubles[1]->name, "bb"));
  inla_cgeneric_vec_tp *bb = data->doubles[1];
  assert(bb->len == nm);

  // prior parameters
  assert(!strcasecmp(data->doubles[2]->name, "prs"));
  inla_cgeneric_vec_tp *prs = data->doubles[2];
  assert(prs->len == 2);

  assert(!strcasecmp(data->doubles[3]->name, "prt"));
  inla_cgeneric_vec_tp *prt = data->doubles[3];
  assert(prt->len == 2);

  assert(!strcasecmp(data->doubles[4]->name, "psigma"));
  inla_cgeneric_vec_tp *psigma = data->doubles[4];
  assert(psigma->len == 2);

  assert(!strcasecmp(data->mats[0]->name, "tt"));
  inla_cgeneric_mat_tp *tt = data->mats[0];
  assert(tt->nrow == nm);
  assert(tt->ncol == 2);

  nth=0;
  if(iszero(prs->doubles[1])){
    ifix[0] = 1;
  } else {
    ifix[0] = 0;
    nth++;
  }
  if(iszero(prt->doubles[1])) {
    ifix[1] = 1;
  } else {
    ifix[1] = 0;
    nth++;
  }
  if(iszero(psigma->doubles[1])) {
    ifix[2] = 1;
  } else {
    ifix[2] = 0;
    nth++;
  }
  assert(nth<4);
  
  if (theta) {
    // interpretable parameters 
    //  theta = log(range_s, range_t, sigma)
    // g_s = sqrt(8 * v_s) / r_s ;
    //   log(g_s) = c1 + log(r_s) ;
    // g_t = r_t * g_s^a_s / sqrt(8(a_t-1/2)) ;
    //   log(g_t) = c2 + log(r_t) + a_s log(g_s) ;
    // g_e = sqrt(c12 / (sigma * g_t * g_s^{2a-d})) ;
    //   log(g_e) = c3 + log(sigma) + log(g_t) + (d-2a)log(g_s) ;
    ith = 0;
    if(ifix[0]==1){
      lg[0] = cc->doubles[0] - log(prs->doubles[0]);
    } else {
      lg[0] = cc->doubles[0] - theta[ith++];
    }
    if(manifold==0) {
      if (ialpha & (((int)alpha)==1L)) {
	for(int k=0; k<50; k++) {
	  c3 += (2 * (float)k) / (exp(2*lg[0]) + (float)((k*(k+1))));
	}
      }
      if (ialpha & (((int)alpha)==2L)) {
	for(int k=0; k<50; k++) {
	  c3 += (2 * (float)k) / pow2(exp(2*lg[0]) + (float)((k*(k+1))));
	}
      }
      if (ialpha & (((int)alpha)==3L)) {
	for(int k=0; k<50; k++) {
	  c3 += (2 * (float)k) / pow3(exp(2*lg[0]) + (float)((k*(k+1))));
	}
      }
      if (ialpha & (((int)alpha)==4L)) {
	for(int k=0; k<50; k++) {
	  c3 += (2 * (float)k) / pow4(exp(2*lg[0]) + (float)((k*(k+1))));
	}
      }
      if(!ialpha) {
	for(int k=0; k<50; k++) {
	  c3 += (2 * (float)k) / pow(exp(2*lg[0]) + (float)((k*(k+1))), alpha);
	}
      }
    }
    if(ifix[1]==1){
      lg[1] = cc->doubles[1] + alphas*lg[0] + log(prt->doubles[0]);
    } else {
      lg[1] = cc->doubles[1] + alphas*lg[0] + theta[ith++];
    }
    if(ifix[2]==1){
      lg[2] = c3 +aaux*lg[0] -0.5*lg[1] - log(psigma->doubles[0]); 
    } else {
      lg[2] = c3 +aaux*lg[0] -0.5*lg[1] - theta[ith++];
    }
    assert(nth == ith);

    k=0;
    for(i=0; i<nm; i++) {
      a1 = lg[0] * tt->x[k++];
      a2 = lg[1] * tt->x[k++];
      params[i] = exp( 2*(lg[2] + a1 + a2) ) * bb->doubles[i];
    }
  } else {
    for(i=0; i<nm; i++) {
      params[i] = NAN;
    }
  }

  
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }
    
  case INLA_CGENERIC_GRAPH:
    {
      k = 2;
      ret = Calloc(k + 2 * M, double);
      ret[0] = N;        	       /* dimension */
      ret[1] = M;		   /* number of (i <= j) */
      for (i = 0; i < M; i++) {
	ret[k++] = ii->ints[i];
      }
      for (i = 0; i < M; i++) {
	ret[k++] = jj->ints[i];
      }
      break;
    }
    
  case INLA_CGENERIC_Q:
    {

      int offset = 2;
      ret = Calloc(offset + M, double);

      assert(!strcasecmp(data->mats[1]->name, "xx"));
      inla_cgeneric_mat_tp *xx = data->mats[1];
      assert(xx->nrow == nm);
      assert(xx->ncol == M);

      int one=1;
      double zerof = 0.0, onef = 1.0;
      char trans = 'N';
      
      dgemv_(&trans, &M, &nm, &onef, &xx->x[0], &M, params, &one, &zerof, &ret[offset], &one, F_ONE);
      
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
      // return c(P, initials)
      // where P is the number of hyperparameters      
      ret = Calloc(nth+1, double);
      ith = 0;
      ret[ith++] = (double)nth;
      if(ifix[0]==0){
	ret[ith++] = 1.0;
      } 
      if(ifix[1]==0) {
	ret[ith++] = 1.0;
      } 
      if(ifix[2]==0) {
	ret[ith++] = 0.0;
      }
      assert(ith==(nth+1));
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
	double lam1 = -log(prs->doubles[1])/prs->doubles[0];
	ret[0] += log(lam1) - theta[ith] - lam1*exp(-theta[ith]) ;
	ith++;
      }
      if(ifix[1]==0) {
	double lam2 = -log(prt->doubles[1])/prt->doubles[0];      
	ret[0] += log(lam2) - 0.5*theta[ith] - lam2*exp(-0.5*theta[ith]) +log(0.5) ;
	ith++;
      }
      if(ifix[2]==0) {
	double lam3 = -log(psigma->doubles[1])/psigma->doubles[0];
	ret[0] += log(lam3) + theta[ith] - lam3*exp(theta[ith]) ;
	ith++;
      }
      assert(ith==nth);
      /*ret[0] = 0.0;
      double z = 0.0;
      double m[3] = {8, 5, 3};
      double s[3] = {3, 3, 3};
      for(j=0; j<nth; j++) {
	z = (theta[j] - m[j])/s[j];
	ret[0] -= 0.5*SQR(z) + 0.9189385332 + log(s[j]);
      }*/
	
      break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}

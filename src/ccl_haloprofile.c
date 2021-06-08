#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include "ccl.h"

typedef struct {
  double label;
  int prof_sign;
  ccl_f1d_t *prof_real;
  double x_min;
  double x_max;
} integ_profiles;


static double profile_projected_integrand(double r_par, void *params)
{
  integ_profiles *p=(integ_profiles *)params;
  double r_perp=p->label;
  double r=sqrt(r_perp*r_perp+r_par*r_par);
  return p->prof_sign*exp(ccl_f1d_t_eval(p->prof_real, log(r)));
}

void ccl_profile_projected(int nr, double *r_arr, double *p_arr, double r_min, double r_max,
			   int nr_perp, double *r_perp_arr, double *p_out)
{
  int ir;
  integ_profiles *ip=(integ_profiles *)malloc(sizeof(integ_profiles));
  // TODO check memory everywhere
  ip->prof_sign = 1;
  if(p_arr[0] < 0)
    ip->prof_sign=-1;
  for(ir=0;ir<nr;ir++) {
    p_arr[ir] = log(ip->prof_sign*p_arr[ir]);
    r_arr[ir] = log(r_arr[ir]);
  }
  int status=0;
  ip->prof_real=ccl_f1d_t_new(nr, r_arr, p_arr, p_arr[0], p_arr[nr-1],
			      ccl_f1d_extrap_linx_liny, ccl_f1d_extrap_linx_liny, &status);
  for(ir=0;ir<nr;ir++) {
    p_arr[ir] = ip->prof_sign*exp(p_arr[ir]);
    r_arr[ir] = exp(r_arr[ir]);
  }
  ip->x_min=r_min;
  ip->x_max=r_max;

  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(100000);
  F.function = &profile_projected_integrand;
  F.params = ip;
  for(ir=0;ir<nr;ir++) {
    int gslstatus;
    double res, error, rmax;
    rmax = fmax(ip->x_max, 100*r_perp_arr[ir]);
    ip->label = r_perp_arr[ir];
    //gslstatus = gsl_integration_qag(&F, ip->x_min, rmax, 0, 1E-5, 10000,
    //                                GSL_INTEG_GAUSS61, w, &res, &error);
    gslstatus = gsl_integration_qagiu(&F, 0., 0., 1E-10, 100000,
				      w, &res, &error);
    if(gslstatus)
      printf("%d %d %lE %lE\n", ir, gslstatus, res, error);
    p_out[ir] = 2*res;
  }
  gsl_integration_workspace_free(w);

  ccl_f1d_t_free(ip->prof_real);
  free(ip);
}

static double einasto_norm_integrand(double x, void *params)
{
  double alpha = *((double *)(params));
  return x*x*exp(-2*(pow(x,alpha)-1)/alpha);
}

void ccl_einasto_norm_integral(int n_m, double *r_s, double *r_delta, double *alpha,
			       double *norm_out,int *status)
{
#pragma omp parallel default(none)			\
  shared(n_m, r_s, r_delta, alpha, norm_out, status)
  {
    int ii;
    int status_this=0;
    gsl_function F;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    
    if (w == NULL)
      status_this = CCL_ERROR_MEMORY;
    
    if(status_this == 0) {
#pragma omp for
      for(ii=0;ii<n_m;ii++) {
	int qagstatus;
	double result, eresult;
	double x_max = r_delta[ii]/r_s[ii];
	F.function = &einasto_norm_integrand;
	F.params = &(alpha[ii]);
	qagstatus = gsl_integration_qag(&F, 0, x_max, 0, 1E-4,
					1000, GSL_INTEG_GAUSS31,
					w, &result, &eresult);
	if(qagstatus != GSL_SUCCESS) {
	  ccl_raise_gsl_warning(qagstatus, "ccl_haloprofile.c: ccl_einasto_norm_integral():");
	  status_this = CCL_ERROR_INTEG;
	  result = NAN;
	}
	norm_out[ii] = 4 * M_PI * r_s[ii] * r_s[ii] * r_s[ii] * result;
      }
    } //end omp for
  
    gsl_integration_workspace_free(w);
    if(status_this) {
      #pragma omp atomic write
      *status = status_this;
    }
  } //end omp parallel
}

static double hernquist_norm_integrand(double x, void *params)
{
  double opx=1+x;
  return x*x/(x*opx*opx*opx);
}

void ccl_hernquist_norm_integral(int n_m, double *r_s, double *r_delta,
			       double *norm_out,int *status)
{
#pragma omp parallel default(none)		\
  shared(n_m, r_s, r_delta, norm_out, status)
  {
    int ii;
    int status_this=0;
    gsl_function F;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    
    if (w == NULL)
      status_this = CCL_ERROR_MEMORY;
    
    if(status_this == 0) {
#pragma omp for
      for(ii=0;ii<n_m;ii++) {
	int qagstatus;
	double result, eresult;
	double x_max = r_delta[ii]/r_s[ii];
	F.function = &hernquist_norm_integrand;
	F.params = NULL;
	qagstatus = gsl_integration_qag(&F, 0, x_max, 0, 1E-4,
					1000, GSL_INTEG_GAUSS31,
					w, &result, &eresult);
	if(qagstatus != GSL_SUCCESS) {
	  ccl_raise_gsl_warning(qagstatus, "ccl_haloprofile.c: ccl_hernquist_norm_integral():");
	  status_this = CCL_ERROR_INTEG;
	  result = NAN;
	}
	norm_out[ii] = 4 * M_PI * r_s[ii] * r_s[ii] * r_s[ii] * result;
      }
    } //end omp for
  
    gsl_integration_workspace_free(w);
    if(status_this) {
      #pragma omp atomic write
      *status = status_this;
    }
  } //end omp parallel
}

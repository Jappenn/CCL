/** @file */

#ifndef __CCL_HALOPROFILE_H_INCLUDED__
#define __CCL_HALOPROFILE_H_INCLUDED__

CCL_BEGIN_DECLS


/**
 * Computes the normalization of the Einasto profile.
 * @param cosmo: cosmology object containing parameters
 * @param n_m: number of elements of the input and output arrays.
 * @param r_s: scale radius
 * @param r_delta: virial radius
 * @param alpha: index
 * @param norm_out: output array
 * @param status: Status flag: 0 if there are no errors, non-zero otherwise
 * @return void
 */
void ccl_einasto_norm_integral(int n_m, double *r_s, double *r_delta, double *alpha,
			       double *norm_out,int *status);


/**
 * Computes the normalization of the Hernquist profile.
 * @param cosmo: cosmology object containing parameters
 * @param n_m: number of elements of the input and output arrays.
 * @param r_s: scale radius
 * @param r_delta: virial radius
 * @param norm_out: output array
 * @param status: Status flag: 0 if there are no errors, non-zero otherwise
 * @return void
 */
void ccl_hernquist_norm_integral(int n_m, double *r_s, double *r_delta,
				 double *norm_out,int *status);

void ccl_profile_projected(int nr, double *r_arr, double *p_arr, double r_min, double r_max,
			   int nr_perp, double *r_perp_arr, double *p_out);

CCL_END_DECLS
#endif

/*
 ============================================================================
 Name        : FXPF.c
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : 
 Description : 	This program determines the Fixed-Point Format of a
 			   	LTI filter given in state-space representation:
 			   		|x(k+1) = Ax(k) + Bu(k)
 			   		|y(k)   = Cx(k) + Du(k)
				where
					*u(k) is an input vector of size q x 1
					*x(k) is internal state vector of size n x 1
					*y(k) is the filter output vector of size p x 1
				Input of the algorithm:
					* ur such that for all k |u(k)| <= ur,
					  ur is a 1 x q vector
					* w_u, w_x, w_y the restrictions on the wordlengths
						of vectors u(k), x(k), y(k) respectively
					  w_u is a 1 x q vector
					  w_x is a n x 1 vector
					  w_y is a p x 1 vector


 ============================================================================
 */

#include "fxpf.h"
#include <time.h>

/*Allocates the space for the result structure fields.
	If length is zero, exit with error.	

 */
void fxpf_result_allocate(fxpf_result *result, uint64_t length )
{
	if(!length)
	{
		fprintf(stderr, "ERROR: Trying to initialize fxpf_result with length = 0.\n");
		return;
	}
	result->length = length;
	result->error = allocateMPFIMatrix(1, length, 64);
	result->w = (uint64_t*)calloc(length, length * sizeof(uint64_t));
	result->msb = (uint64_t*)calloc(length, length * sizeof(uint64_t));
	result->lsb = (int*)calloc(length, length * sizeof(int));
}

/* Frees the space by the pointer of the fxpf_result. */
void fxpf_result_free(fxpf_result *result)
{
	if(!result->w) free(result->w);
	if(!result->msb) free(result->msb);
	if(!result->lsb) free(result->lsb);
	freeMPFIMatrix(result->error, 1, result->length);
}

/* Copies the values of an array of wordlength constraints
w (which muct be non-null pointer for an array of length elements)
to the fxpf_result structure. 
If w is a null pointer, print an error.

 */
void fxpf_result_setW(fxpf_result *result, uint64_t *w)
{
	if(!w) 
	{
		fprintf(stderr, "ERROR: Null pointer occured while trying to initialize wordlength constraints.\n");
		return;
	}
	int i;
	for(i = 0; i < result->length; ++i)
	{
		if(w[i] < 4 || isnan(w[i]))
		{
			fprintf(stderr, "ERROR: wordlength constraint must be >=4 bits but it is w[%i] = %llu \n", i, w[i]);
			return;
		}
		result->w[i] = w[i];
	}

}

/* Copies the values of an array of MSB msb (which muct be non-null pointer for an array of length elements)
to the fxpf_result structure. 
If msb is a null pointer, print an error.
*/

void fxpf_result_setMSB(fxpf_result *result, uint64_t *msb)
{
	if(!msb) 
	{
		fprintf(stderr, "ERROR: Null pointer occured while trying to copy MSBs to fxpf_result.\n");
		return;
	}
	UINT64MatrixCopy(result->msb, msb, 1, result->length);
}

/* Copies the values of an array of LSBs lsb (which muct be non-null pointer for an array of length elements)
to the fxpf_result structure. 
If lsb is a null pointer, print an error.
*/
void fxpf_result_setLSB(fxpf_result *result, int *lsb)
{
	if(!lsb) 
	{
		fprintf(stderr, "ERROR: Null pointer occured while trying to copy LSBs to fxpf_result.\n");
		return;
	}
	INTMatrixCopy(result->lsb, lsb, 1, result->length);
}

/* Copies the values of an array of error-vector (which muct be non-null pointer for an array of length elements)
to the fxpf_result structure. 
If error is a null pointer, print an error.
*/
void fxpf_result_setError(fxpf_result *result, mpfi_t *error)
{
	if(!error) 
	{
		fprintf(stderr, "ERROR: Null pointer occured while trying to copy error vector to fxpf_result.\n");
		return;
	}
	int i;
	mp_prec_t prec;
	for(i = 0; i < result->length; ++i)
	{
		prec = mpfi_get_prec(error[i]) > mpfi_get_prec(result->error[i]) ? mpfi_get_prec(error[i]) : mpfr_get_default_prec();
		mpfi_set_prec(result->error[i], prec);

		mpfi_set(result->error[i], error[i]);
	}
}

/* Sets the number of additional steps. The integer value steps
must be positive or zero, otherwise print an error.
 */
void fxpf_result_setAdditionalSteps(fxpf_result *result, int steps)
{
	if(steps < 0) fprintf(stderr, "WARNING: Trying to set negative numver of additional steps. \n");
	result->additionalSteps = steps;
}

/* Prints results in format
	wordlengths
	msb
	lsb
	errors
	additional steps
	(empty line)
 */
void fxpf_result_print(FILE *file, fxpf_result *result)
{
	// fprintf(file, "Printing result:\n");
	// fprintf(file, "%llu\n",result->length);
	UINT64MatrixPrint(file, result->w, 1, result->length);
	// fprintf(file, "\n");
	// UINT64MatrixPrint(file, result->msb, 1, result->length);
	// fprintf(file, "\n");
	// INT64MatrixPrint(file, result->lsb, 1, result->length);
	int i;
	double d;

	for(i = 0; i < result->length; ++i)
	{
		d = mpfi_get_d(result->error[i]);
		fprintf(file, "%e \t", d);

	}
	fprintf(file, "\n");	
}



/* Computes the precision of WCPG computation of filter H such that
the MSB computed with above mentioned WCPG is either exact, either overestimated by one.
The precision eps is computed via following formula:
	eps = max_i {eps_i}
	eps_i <= (WCPGLowerBound(H) * ur)_i / alpha * sum_j {ur_j} for i=1:p

with alpha = 1 or 2.
Note: we are going to compute WCPG lower bound with precision 64.

Returns non-zero if epsilon was computied; and zero otherwise.
 */
int getWCPGeps(mpfr_t epsWCPG, filter *H, uint64_t alpha)
{

	mp_exp_t res;
	mp_prec_t prec = 64;

	/* computing the WCPG lower bound matrix */
	mp_exp_t epsWCPGLowerBound = -prec;

	mpfr_t *WLowerBound;
	WLowerBound = allocateMPFRMatrix(H->p, H->q, prec);

	// printf("Filter: \n");
	// filter_print(stderr, H);
	if(!WCPGLowerBound_double(WLowerBound, H->A, H->B, H->C, H->D, H->n, H->p, H->q, -64))
	{
		fprintf(stderr, "Could not compute lower bound on the Worst Case Peak Gain matrix. \n");
		freeMPFRMatrix(WLowerBound, H->p, H->q);
		return 0;
	}

	mpfi_t *eps;
	eps = allocateMPFIMatrix(H->p, 1, prec);

	mpfi_t tmp;
	mpfi_init(tmp);

	/* computing denumerator sum_j (ur_j) */
	mpfi_t sum_ur;
	mpfi_init_set_ui(sum_ur, 0);
	int i, j;
	for(i = 0; i < H->q; ++i)
	{
		mpfi_add_d(sum_ur, sum_ur, H->ur[i]);
	}

	 // fprintf(stderr, "Denum sum_j U_r: \n");
	 // mpfi_out_str(stderr, 10, 10, sum_ur);

	/* computing the epsilon vector */
	mpfi_t *WCPG_ur;
	WCPG_ur = allocateMPFIMatrix(H->p, 1, prec);
	for(i = 0; i < H->p; ++i)
	{
		mpfi_set_ui(WCPG_ur[i], 0);
		for(j = 0; j < H->q; ++j)
		{

			mpfi_set_fr(tmp, WLowerBound[i * H->q + j]);

			mpfi_mul_d(tmp, tmp, H->ur[j]);
	
			mpfi_add(WCPG_ur[i], WCPG_ur[i], tmp);
			
		}
		mpfi_div(eps[i], WCPG_ur[i],sum_ur); 
		if(mpfi_nan_p(eps[i]))
		{
			fprintf(stderr, "\n Could not compute the epsWCPG: eps[%d] is NaN \n", i);

			freeMPFRMatrix(WLowerBound, H->p, H->q);
			mpfi_clear(tmp);
			mpfi_clear(sum_ur);
			freeMPFIMatrix(eps, 1, H->p);
			freeMPFIMatrix(WCPG_ur, H->p, 1);
			return 0;
		}
		 // printf("\n eps[%d] = \n", i);
		 // mpfi_out_str(stderr, 10, 5, eps[i]);
	}
	

	mpfr_t curr, min;
	mpfr_init_set_ui(curr, 1, MPFR_RNDN);
	mpfr_init(min);							//setting min
	mpfi_get_left(min, eps[0]);

	int successflag = 1;

	if(alpha == (uint64_t)2)
	{
		for(i = 0; i < H->p; ++i)
		{
			mpfi_div_ui(eps[i], eps[i], 2);
			mpfi_get_left(curr, eps[i]);
			if(mpfr_cmp(min, curr)) mpfr_set(min, curr, MPFR_RNDD);
		}
		mpfr_div_ui(epsWCPG, min, 4, MPFR_RNDN);	//just adding additional 2 bits
		
		if(mpfr_nan_p(epsWCPG))
		{
			fprintf(stderr, "\n epsWCPG is NaN \n");
			successflag = 0;
		}

	}
	if (alpha == (uint64_t)1)
	{
		for(i = 0; i < H->p; ++i)
		{
			mpfi_get_left(curr, eps[i]);
			if(mpfr_cmp(min, curr)) mpfr_set(min, curr, MPFR_RNDD);
		}
		mpfr_div_ui(epsWCPG, min, 4, MPFR_RNDN);	//just adding additional 2 bits
		
		if(mpfr_nan_p(epsWCPG))
		{
			fprintf(stderr, "\n epsWCPG is NaN \n");
			successflag = 0;
		}

	}


	mpfr_clear(curr);
	mpfr_clear(min);
	freeMPFRMatrix(WLowerBound, H->p, H->q);
	mpfi_clear(tmp);
	mpfi_clear(sum_ur);
	freeMPFIMatrix(eps, 1, H->p);
	freeMPFIMatrix(WCPG_ur, H->p, 1);

	return successflag;
}

/*  Computes MSB for the Step 1 of the algorithm: a modification Hz1 of filter is 
	considered instead of filter H. This filter Hz1 does not take into account 
	any computational errors.
	Given wordlength constraints, the MSBs of this filter are evaluated with
		m_i = ceil[log2(WCPG(Hz1)*u_bound) - log2(A - 2^(1-w_i)) ]	(1)

	If the WCPG of filter Hz1 is computed with epsWCPG satisfying condifiton
		epsWCPG   <= min_i {WCPGLowerBound(Hz1)*u_bound / 2*sum_j u_bound_j }

	then, the difference between the exact MSBs and the MSBs computed with (1)
	is 0 or 1.

	Inputs:
		* filter Hz1 
		* wordlengths w
	Outputs:
		* msb_implemented

	Returns non-zero value in case of success and zero otherwise
   */
int msb_step1(mpfi_t *msb, uint64_t *msb_left, uint64_t *msb_right, filter *H, uint64_t *w)
{
	
	mp_prec_t defaultPrec = 64;

	mpfr_t epsWCPG;
	mpfr_init2(epsWCPG, defaultPrec);

	if(!getWCPGeps(epsWCPG, H, (uint64_t)1) || mpfr_cmp_ui(epsWCPG, (uint64_t)0) <= 0 )
	{
		fprintf(stderr, "Could not determine the epsilon for WCPG computation. Exit with error. \n");
		mpfr_clear(epsWCPG);
		return 0;
	}

	// printf("\n ...epsWCPG: \n");
	// mpfr_out_str(stderr, 10, 20, epsWCPG, MPFR_RNDN);

	mpfr_t *WCPG;
	WCPG = allocateMPFRMatrix(H->p, H->q, defaultPrec);
	
	if (!WCPG_ABCD_mprec(WCPG, H->A, H->B, H->C, H->D, H->n, H->p, H->q, epsWCPG))
	{
		fprintf(stderr, "Could not compute WCPG. Exit with error\n");
		freeMPFRMatrix(WCPG, H->p, H->q);
		mpfr_clear(epsWCPG);
		return 0;
	}	
	// fprintf(stderr, "\n WCPG matrix: \n");
	// writeMPFRMatrix(stderr, WCPG, H->p, H->q, 10, MPFR_RNDN);
	

	mpfr_t tmp_mpfr;
	mpfr_init(tmp_mpfr);
	mp_prec_t WCPGprec = getMaxPrecMPFRMatrix(WCPG, H->p, H->q);
	defaultPrec = WCPGprec >= defaultPrec ? WCPGprec : defaultPrec;		//set default precision to the max precision of the WCPG matrix or leave it unchanged 


	mpfi_t WCPG_ur, tmp;
	mpfi_init2(WCPG_ur, 2 * defaultPrec); 	//we will store [(<<H>> * ur)_i] in this variable
	mpfi_set_ui(WCPG_ur, 0);
	mpfi_init2(tmp, 2 * defaultPrec); 
	mpfi_set_ui(tmp, 0);

	mpfi_t WCPG_ij_interval;					//we will store interval values of the WCPG_ij
	mpfi_init2(WCPG_ij_interval, defaultPrec);	//with precision being at least the max precision of the WCPG matrix

	mpfi_t log2WCPG_ur, log2const;
	mpfi_init2(log2WCPG_ur, defaultPrec);
	mpfi_init2(log2const, defaultPrec);
	int i, j;
	for(i = 0; i < H->p; ++i)
	{
		/* Compute [log2[(<<H>> * ur)_i]] */
		mpfi_set_ui(WCPG_ur, 0);
		for(j = 0; j < H->q; ++j)
		{
			mpfi_set_fr(WCPG_ij_interval, WCPG[i * H->q + j]);		// convert WCPG_ij to interval
			/* Here precision of H->ur[j] is at most defaultprec.
			Setting precision of tmp equal to 2 * defaultprec guarantees the 
			interval bounds to be exact.  */
			mpfi_mul_d(tmp, WCPG_ij_interval, H->ur[j]);			
			mpfi_add(WCPG_ur, WCPG_ur, tmp);						//WCPG_ur = [(<<H>> * ur)_i]
		}
		mpfi_log2(log2WCPG_ur,WCPG_ur);
		if(mpfi_nan_p(log2WCPG_ur))
		{
			fprintf(stderr, "\n Could not compute MSB: log2WCPG_ur is NaN \n");
			mpfi_clear(WCPG_ur);
			mpfi_clear(tmp);
			mpfi_clear(log2WCPG_ur);
			mpfi_clear(log2const);
			freeMPFRMatrix(WCPG, H->p, H->q);
			mpfr_clear(tmp_mpfr);
			mpfr_clear(epsWCPG);
			mpfi_clear(WCPG_ij_interval);
			return 0;

		}
		  // fprintf(stderr, "\nStep 1:log2WCPG_ur_%d: \n", i);
		  // mpfi_out_str(stderr, 10, 10, log2WCPG_ur);

		/* Compute log2const = [log2(1 - 2^(1 - w_j))] */
		mpfi_set_prec(tmp, defaultPrec);
		mpfi_set_si(tmp, 1 - w[i]);
		if(mpfi_cmp_ui(tmp, 0) == 0)	
		{
			/* Here w[i] is <= 1, which is violation of input constraints. */
			fprintf(stderr, "ERROR: wordlength is 1: 1 - w = %f.\n", mpfi_get_d(tmp));
			mpfi_clear(WCPG_ur);
			mpfi_clear(tmp);
			mpfi_clear(log2WCPG_ur);
			mpfi_clear(log2const);
			freeMPFRMatrix(WCPG, H->p, H->q);
			mpfr_clear(tmp_mpfr);
			mpfr_clear(epsWCPG);
			mpfi_clear(WCPG_ij_interval);
		}

		mpfi_exp2(tmp, tmp);
		mpfi_neg(tmp,tmp);
		mpfi_log1p(log2const, tmp);

		if(mpfi_nan_p(log2const))
		{
			fprintf(stderr, "\n ERROR: Could not compute MSB: log2(1 - 2^(1-w)) is NaN \n");
			mpfi_clear(WCPG_ur);
			mpfi_clear(tmp);
			mpfi_clear(log2WCPG_ur);
			mpfi_clear(log2const);
			freeMPFRMatrix(WCPG, H->p, H->q);
			mpfr_clear(tmp_mpfr);
			mpfr_clear(epsWCPG);
			mpfi_clear(WCPG_ij_interval);
			return 0;
		}
		 // fprintf(stderr, "\n log2const: \n");
		 // mpfi_out_str(stderr, 10, 10, log2const);

		/* Compute [log2[(<<H>> * ur)_i]] - [log2(1 - 2^(1 - w_j))] */
	
		mpfi_set_prec(msb[i], mpfi_get_prec(log2WCPG_ur));
		mpfi_sub(msb[i], log2WCPG_ur, log2const);
		
		if(mpfi_nan_p(msb[i]))
		{
			fprintf(stderr, "\n ERROR: Could not compute MSB: msb[%d] is NaN \n", i);
			mpfi_clear(WCPG_ur);
			mpfi_clear(tmp);
			mpfi_clear(log2WCPG_ur);
			mpfi_clear(log2const);
			freeMPFRMatrix(WCPG, H->p, H->q);
			mpfr_clear(tmp_mpfr);
			mpfr_clear(epsWCPG);
			mpfi_clear(WCPG_ij_interval);
			return 0;
		}

		mpfr_set_prec(tmp_mpfr, 2 * defaultPrec);

		if(mpfi_get_left(tmp_mpfr, msb[i]) < 0)
		{
			fprintf(stderr, "Could not correctly round inf of the MSB interval. \n");
		}

		if(mpfr_ceil(tmp_mpfr, tmp_mpfr) < 0)
		{
			fprintf(stderr, "Could not correctly round ceil(inf(MSB)). \n");
		}
		msb_left[i] = mpfr_get_ui(tmp_mpfr, MPFR_RNDD);


		mpfr_set_prec(tmp_mpfr, 2 * defaultPrec);
		if(mpfi_get_right(tmp_mpfr, msb[i]) < 0)
		{
			fprintf(stderr, "Could not correctly round sup of the MSB interval. \n");
		}
		if(mpfr_ceil(tmp_mpfr,tmp_mpfr) < 0)
		{
			fprintf(stderr, "Could not correctly round ceil(sup(MSB)). \n");
		}
		msb_right[i] = mpfr_get_ui(tmp_mpfr, MPFR_RNDU);
		
	}

		freeMPFRMatrix(WCPG, H->p, H->q);
		mpfr_clear(epsWCPG);
		mpfi_clear(WCPG_ur);
		mpfi_clear(tmp);
		mpfr_clear(tmp_mpfr);
		mpfi_clear(log2WCPG_ur);
		mpfi_clear(log2const);
		mpfi_clear(WCPG_ij_interval);

return 1;

}


/*  Computes MSB for the Step 2 of the algorithm: output of the 
	implemented filter Hz1 is represented as sum of the outputs of the
	filter Hz1, which does not take into account computational errors,
	and error-filter Delta.
	Then, MSB of the implemented filter, taking into account the computational
	errors is evaluated with:

		m_i = ceil[log2{WCPG(Hz1)*u_bound)_i + WCPG(Delta)*eps_bound)_i } - log2{A - 2^(1-w_i)} ]   (1)

	If the WCPGs of filters are computes with epsWCPG satisfying condifiton
		epsWCPG_Hz1   <= min_i {WCPGLowerBound(Hz1)*u_bound / 2*sum_j u_bound_j }
		epsWCPG_Delta <= min_i {WCPGLowerBound(Delta)*u_bound / 2*sum_j u_bound_j }

	then, the difference between the exact MSBs and the MSBs computed with (1)
	is 0 or 1.

	Inputs:
		* filter Hz1 
		* filter Delta
		* msb_Hz1
	Outputs:
		* msb_implemented
   */

int msb_step2(mpfi_t *msb, uint64_t *msb_left, uint64_t *msb_right, filter *Hz1, filter *Delta, uint64_t *w, uint64_t *lsb_max, fxpf_result *result)
{
	mp_prec_t defaultPrec = 64;

	mpfr_t epsWCPG_Hz1;
	mpfr_init2(epsWCPG_Hz1, defaultPrec);

	if(!getWCPGeps(epsWCPG_Hz1, Hz1, (uint64_t)2) || mpfr_cmp_ui(epsWCPG_Hz1, (uint64_t)0) <= 0 )
	{
		fprintf(stderr, "Could not determine the epsilon for WCPG_Hz1 computation. Exit with error. \n");
		mpfr_clear(epsWCPG_Hz1);
		return 0;
	}
	

	mpfr_t epsWCPG_Delta;
	mpfr_init2(epsWCPG_Delta, defaultPrec);

	if(!getWCPGeps(epsWCPG_Delta, Delta, (uint64_t)2) || mpfr_cmp_ui(epsWCPG_Delta, (uint64_t)0) <= 0 )
	{
		fprintf(stderr, "Could not determine the epsilon for WCPG_Delta computation. Exit with error. \n");
		mpfr_clear(epsWCPG_Hz1);
		mpfr_clear(epsWCPG_Delta);
		return 0;
	}

	 // printf("\n ...epsWCPG_Hz1: \n");
	 // mpfr_out_str(stderr, 10, 20, epsWCPG_Hz1, MPFR_RNDN);
	 // printf("\n ...epsWCPG_Delta: \n");
	 // mpfr_out_str(stderr, 10, 20, epsWCPG_Delta, MPFR_RNDN);

	/* Computing WCPGs */
	
	mpfr_t *WCPG_Hz1;
	WCPG_Hz1 = allocateMPFRMatrix(Hz1->p, Hz1->q, defaultPrec);

	if (!WCPG_ABCD_mprec(WCPG_Hz1, Hz1->A, Hz1->B, Hz1->C, Hz1->D, Hz1->n, Hz1->p, Hz1->q, epsWCPG_Hz1))
	{
		fprintf(stderr,"\nERROR: Could not compute WCPG_Hz1. Exiting.\n");
		freeMPFRMatrix(WCPG_Hz1,Hz1->p, Hz1->q);
		mpfr_clear(epsWCPG_Hz1);
		mpfr_clear(epsWCPG_Delta);
		return 0;
	}	


 	 
	mpfr_t *WCPG_Delta;
	WCPG_Delta = allocateMPFRMatrix(Delta->p, Delta->q, defaultPrec);
	if (!WCPG_ABCD_mprec(WCPG_Delta, Delta->A, Delta->B, Delta->C, Delta->D, Delta->n, Delta->p, Delta->q, epsWCPG_Delta))
	{
		fprintf(stderr, "\nERROR: Could not compute WCPG_Delta. Exiting.\n");
		freeMPFRMatrix(WCPG_Delta,Delta->p, Delta->q);
		freeMPFRMatrix(WCPG_Hz1,Hz1->p, Hz1->q);
		mpfr_clear(epsWCPG_Hz1);
		mpfr_clear(epsWCPG_Delta);
		return 0;
	}	
	
	 //fprintf(stderr, "\n <---------Debug-------->\nWCPG_Hz1 matrix: \n");
	 //writeMPFRMatrix(stderr, WCPG_Hz1, Hz1->p, Hz1->q, 10, MPFR_RNDN);

	 //fprintf(stderr, "\n <---------Debug-------->\nWCPG_Delta matrix: \n");
	 //writeMPFRMatrix(stderr, WCPG_Delta, Delta->p, Delta->q, 10, MPFR_RNDN);

	
	

	
	mp_prec_t precDelta = getMaxPrecMPFRMatrix(WCPG_Delta, Delta->p, Delta->q);
	mp_prec_t precHz1 = getMaxPrecMPFRMatrix(WCPG_Hz1, Hz1->p, Hz1->q);
	defaultPrec = precHz1 >= precDelta ? precHz1 : precDelta;

	mpfi_t tmp_mpfr;
	mpfi_init2(tmp_mpfr, defaultPrec);

	/* Computing element by element the MSB:
		m_i = ceil[log2{WCPG(Hz1)*u_bound)_i + WCPG(Delta)*eps_bound)_i } - log2{A - 2^(1-w_i)} ]
	*/


	mpfi_t WCPG_ij_interval;					//we will store interval values of the WCPG_ij
	mpfi_init2(WCPG_ij_interval, defaultPrec);	//with precision being at least the max precision of the WCPG matrices

	mpfi_t WCPG_Hz1_ur, tmp;
	mpfi_init2(WCPG_Hz1_ur, 2 * defaultPrec);
	mpfi_set_ui(WCPG_Hz1_ur, 0);

	mpfi_init2(tmp, 2 * defaultPrec); 
	mpfi_set_ui(tmp, 0);

	mpfi_t WCPG_Delta_eps;
	mpfi_init2(WCPG_Delta_eps, 2 * defaultPrec); 
	mpfi_set_ui(WCPG_Delta_eps, 0);
	
	mpfi_t log2WCPGs, log2const;
	mpfi_init2(log2WCPGs, 2 * defaultPrec);
	mpfi_init2(log2const, defaultPrec);
	int i, j;

	mpfi_t *Deltaeps_vector;
	Deltaeps_vector = allocateMPFIMatrix(Hz1->p, 1, 2 * defaultPrec);
	mpfi_t *Hz1ur_vector;
	Hz1ur_vector = allocateMPFIMatrix(Hz1->p, 1, 2 * defaultPrec);

	/* Compute the interval vectors 
			Hz1ur_vector = [<<Hz1>> * u ] 
			Deltaeps_vector = [(<<Delta>> * eps)]
	*/
	for(i = 0; i < Hz1->p; ++i)
	{
		mpfi_set_ui(WCPG_Hz1_ur, 0);
		mpfi_set_ui(WCPG_Delta_eps, 0);

		/* Compute [(<<Hz1>> * ur)_i] = sum_{j=1}^{q} [<<Hz1>>_ij * u_j]*/
		for(j = 0; j < Hz1->q; ++j)
		{
			/* convert WCPG_ij to interval */
			mpfi_set_fr(WCPG_ij_interval, WCPG_Hz1[i * Hz1->q + j]);		
			// printf("\n ...Multiplying the WCPG_Hz1[%d] = \n", i * Hz1->q + j);
			// mpfi_out_str(stderr, 10, 10,WCPG_ij_interval);
			// printf("\t by ur[%d] = %d \n", j,Hz1->ur[j]);
			mpfi_set_ui(tmp, 0);

			/* compute the interval [tmp] = <<H_z>>_ij * u_j */
			mpfi_mul_d(tmp, WCPG_ij_interval, Hz1->ur[j]);

			/* add the [tmp] interval to the [(<<Hz1>> * ur)_i] */
			mpfi_add(WCPG_Hz1_ur, WCPG_Hz1_ur, tmp);				
		}

		/* copy the i-th element of the result vector to variable Hz1ur_vector */
		mpfi_set(Hz1ur_vector[i], WCPG_Hz1_ur);
		

		// printf("\n ++Step i = %d. Hz1ur_vector[i]: \n", i);
		// mpfi_out_str(stderr, 10, 10, Hz1ur_vector[i]);

		
		/* Compute [(<<Delta>> * eps)_i] = sum_{j=1}^{q} [<<Delta>>_ij * eps_j]*/
		for(j = 0; j < Delta->q; ++j)
		{
			/* convert WCPG_ij to interval */
			mpfi_set_fr(WCPG_ij_interval, WCPG_Delta[i * Delta->q + j]);		
			// printf("\nIteration %d: WCPG_Delta[%d * %d + %d] = ",i,i,Delta->q, j);
			// mpfi_out_str(stderr, 10, 10, tmp); 
			mpfi_set_ui(tmp, 0);
			/* compute the interval [tmp] = <<Delta>>_ij * eps_j */
			mpfi_mul_d(tmp, WCPG_ij_interval, Delta->ur[j]);
			// printf("\nIteration %d: WCPG_Delta[%d * %d + %d] * %f = \n",i,i,Delta->q, j, j, Delta->ur[j]);
			// mpfi_out_str(stderr, 10, 10, tmp);

			/* add the [tmp] interval to the [(<<Delta>> * eps)_i] */
			mpfi_add(WCPG_Delta_eps, WCPG_Delta_eps, tmp);					
		}
		/* copy the i-th element of the result vector to variable Deltaeps_vector */
		mpfi_set(Deltaeps_vector[i], WCPG_Delta_eps);
		
		// printf("\n --Step i = %d. Deltaeps_vector[i]: \n", i);
		// mpfi_out_str(stderr, 10, 10, WCPG_Delta_eps);

		/* Compute [<<Hz1>> * u ] + [(<<Delta>> * eps)] */
		mpfi_add(log2WCPGs, WCPG_Hz1_ur, WCPG_Delta_eps);

	
		/* Compute log2{(<<Hz1>> * ur)_i]] + [(<<Delta>> * eps)_i]} */
		mpfi_log2(log2WCPGs, log2WCPGs);

		if(mpfi_nan_p(log2WCPGs))
		{
			fprintf(stderr, "\n Could not compute MSB: log2WCPGs is NaN \n");
			mpfr_clear(tmp_mpfr);
			mpfi_clear(WCPG_Hz1_ur);
			mpfi_clear(tmp);
			mpfi_clear(WCPG_Delta_eps);
			mpfi_clear(log2WCPGs);
			mpfi_clear(log2const);
			freeMPFIMatrix(Deltaeps_vector, Hz1->p, 1);
			freeMPFIMatrix(Hz1ur_vector, Hz1->p, 1);
			freeMPFRMatrix(WCPG_Delta,Delta->p, Delta->q);
			freeMPFRMatrix(WCPG_Hz1,Hz1->p, Hz1->q);
			mpfr_clear(epsWCPG_Hz1);
			mpfr_clear(epsWCPG_Delta);
			mpfi_clear(WCPG_ij_interval);
			return 0;
		}


		/* Compute log2const = [log2(1 - 2^(1 - w_j))] */
		mpfi_set_prec(tmp, defaultPrec);
		mpfi_set_si(tmp, 1 - w[i]);
		mpfi_exp2(tmp, tmp);
		mpfi_neg(tmp,tmp);
		mpfi_log1p(log2const, tmp);

		if(mpfi_nan_p(log2const))
		{
			fprintf(stderr, "\n Could not compute MSB: log2const is NaN \n");
			mpfr_clear(tmp_mpfr);
			mpfi_clear(WCPG_Hz1_ur);
			mpfi_clear(tmp);
			mpfi_clear(WCPG_Delta_eps);
			mpfi_clear(log2WCPGs);
			mpfi_clear(log2const);
			freeMPFIMatrix(Deltaeps_vector, Hz1->p, 1);
			freeMPFIMatrix(Hz1ur_vector, Hz1->p, 1);
			freeMPFRMatrix(WCPG_Delta,Delta->p, Delta->q);
			freeMPFRMatrix(WCPG_Hz1,Hz1->p, Hz1->q);
			mpfr_clear(epsWCPG_Hz1);
			mpfr_clear(epsWCPG_Delta);
			mpfi_clear(WCPG_ij_interval);
			return 0;
		}


		/* Compute [log2[(<<Hz1>> * ur)_i]] + [(<<Delta>> * eps)_i]] - [log2(1 - 2^(1 - w_j))] */
		mpfi_sub(msb[i], log2WCPGs, log2const);

		if(mpfi_nan_p(msb[i]))
		{
			fprintf(stderr, "\n Could not compute MSB: msb[%d] is NaN \n", i);
			mpfr_clear(tmp_mpfr);
			mpfi_clear(WCPG_Hz1_ur);
			mpfi_clear(tmp);
			mpfi_clear(WCPG_Delta_eps);
			mpfi_clear(log2WCPGs);
			mpfi_clear(log2const);
			freeMPFIMatrix(Deltaeps_vector, Hz1->p, 1);
			freeMPFIMatrix(Hz1ur_vector, Hz1->p, 1);
			freeMPFRMatrix(WCPG_Delta,Delta->p, Delta->q);
			freeMPFRMatrix(WCPG_Hz1,Hz1->p, Hz1->q);
			mpfr_clear(epsWCPG_Hz1);
			mpfr_clear(epsWCPG_Delta);
			mpfi_clear(WCPG_ij_interval);
			return 0;
		}
		
		/* Now we compute the mpfr left and right bounds of the MSB interval. */
		mpfr_set_prec(tmp_mpfr, 2 * defaultPrec);

		if(mpfi_get_left(tmp_mpfr, msb[i]) < 0)
		{
			fprintf(stderr, "Could not correctly round inf of the MSB interval. \n");
		}

		if(mpfr_ceil(tmp_mpfr, tmp_mpfr) < 0)
		{
			fprintf(stderr, "Could not correctly round ceil(inf(MSB)). \n");
		}
		msb_left[i] = mpfr_get_ui(tmp_mpfr, MPFR_RNDD);


		mpfr_set_prec(tmp_mpfr, 2 * defaultPrec);

		if(mpfi_get_right(tmp_mpfr, msb[i]) < 0)
		{
			fprintf(stderr, "Could not correctly round sup of the MSB interval. \n");
		}
		if(mpfr_ceil(tmp_mpfr,tmp_mpfr) < 0)
		{
			fprintf(stderr, "Could not correctly round ceil(sup(MSB)). \n");
		}
		msb_right[i] = mpfr_get_ui(tmp_mpfr, MPFR_RNDU);



		/* Check the stopping criteria for the algorithm.
			If the lsb[i] on this iteration of the Step 2 
			has moved more to the right than the lsb_max constraint
			(which is the msb[i] position deduced from the Step 1
			of the algorithm), then the algorithm stops.
		 */
		int64_t lsb_current = msb_right[i] - w[i] + 1;
		int64_t lsb_max_tmp = (int64_t)lsb_max[i];
		if(lsb_current > lsb_max_tmp)
		{
			fprintf(stderr, "Error: LSB reached point when it is larger than initial MSB estimation on Step 1!\n");
			fprintf(stderr, "msb_step1 = %llu and lsb_current = %lld \n",lsb_max[i], lsb_current);
			// fprintf(stderr, "Error!!! LSB reached point when it is larger than wordlength.\n");
			mpfr_clear(tmp_mpfr);
			mpfi_clear(WCPG_Hz1_ur);
			mpfi_clear(tmp);
			mpfi_clear(WCPG_Delta_eps);
			mpfi_clear(log2WCPGs);
			mpfi_clear(log2const);
			freeMPFIMatrix(Deltaeps_vector, Hz1->p, 1);
			freeMPFIMatrix(Hz1ur_vector, Hz1->p, 1);
			freeMPFRMatrix(WCPG_Delta,Delta->p, Delta->q);
			freeMPFRMatrix(WCPG_Hz1,Hz1->p, Hz1->q);
			mpfr_clear(epsWCPG_Hz1);
			mpfr_clear(epsWCPG_Delta);	
			return 0;
		}
		
	}

	//fprintf(stderr, "\n <+++++++++++++++++++++++Debug++++++++++++++++++++++++>\n The bound on the output (eps * WCPG_Delta) [p-1]= \n");
	//mpfi_out_str(stderr, 10, 10, Deltaeps_vector[Hz1->p-1]);
	//fprintf(stderr, "\n <++++++++++++++++++++++++++++++++++++++++++++++++++++>\n");
	

	// fprintf(stderr, "\n The <<Delta>>*eps vector is: \n");
	// writeMPFIMatrix(stderr, Deltaeps_vector, 1, Hz1->p, 10);
	// fprintf(stderr, "\n The <<Hz1>>*ur vector is: \n");
	// writeMPFIMatrix(stderr, Hz1ur_vector, 1, Hz1->p, 10);

	/* Set the error vector in the fxpf_result structure to
	error = <<Delta>> * eps 
	*/
	fxpf_result_setError(result, Deltaeps_vector);

	mpfr_clear(tmp_mpfr);
	mpfi_clear(WCPG_Hz1_ur);
	mpfi_clear(tmp);
	mpfi_clear(WCPG_Delta_eps);
	mpfi_clear(log2WCPGs);
	mpfi_clear(log2const);
	freeMPFIMatrix(Deltaeps_vector, Hz1->p, 1);
	freeMPFIMatrix(Hz1ur_vector, Hz1->p, 1);
	freeMPFRMatrix(WCPG_Delta,Delta->p, Delta->q);
	freeMPFRMatrix(WCPG_Hz1,Hz1->p, Hz1->q);
	mpfr_clear(epsWCPG_Hz1);
	mpfr_clear(epsWCPG_Delta);	
	mpfi_clear(WCPG_ij_interval);	

	return 1;


}


/* Computes the LSB vector given wordlegths and MSB positions in the input.
	Space for lsb vector is assumed preallocated output the function.

	Returns non-zero value if succeed, zero otherwise.*/
int computeLSB(int *lsb, uint64_t *msb, uint64_t *w, uint64_t length)
{
	if(!msb || !w || length < 1)
	{
		fprintf(stderr, "ERROR: Null pointer occured in LSB computation function. \n");
		return 0;
	}
	int i;
	for(i = 0; i < length; ++i)
	{
		lsb[i] = msb[i] - w[i] + 1;
	 fprintf(stderr, "compute LSB: lsb[%d] = msb[i] - w[i] + 1 = %lld - %lld + 1 \n", i, lsb[i], msb[i], w[i] );
	}
	return 1;
}


/* For a LTI IIR filter in state-space representation H={A,B,C,D}, bound on input values 
u and wordlegths for state- and output-variable representation the function, if possible,
determines Fixed-Point Formats for its state- and output-variables such that no overflow 
occurs. <<<<>
Algorithm consists in a loop around step1 and step2 parts.

TODO: finish description

Returns a negative number if failed, otherwise number of additional algorithm 
iterations required to determine FXPF.
 */
int determineFXPF(uint64_t *msb, int *lsb, filter *H, uint64_t *wl, fxpf_result *result)
{
	clock_t begin, end;
	begin = clock();
	// fprintf(stderr, "------------------------------ Computing MSB on step 1 --------------------- \n");
	
	/* Construct modification Hz1 of initial filter H. */
	filter Hz1;
	filter_allocate(&Hz1, H->n, H->p + H->n, H->q);
	createFilterHz1(&Hz1, H);

	// printf("++++++++++++++++Filter H:\n");
	// filter_print(stderr, H);
	// printf("++++++++++++++++Filter Hz1:\n");
	// filter_print(stderr, &Hz1);

	mpfi_t *msb_Hz1;
	msb_Hz1 = allocateMPFIMatrix(Hz1.p, 1, 64);
	uint64_t *msb_inf_Hz1;
	msb_inf_Hz1 = (uint64_t*)calloc(Hz1.p, Hz1.p * sizeof(uint64_t)); 
	uint64_t *msb_sup_Hz1;
	msb_sup_Hz1 = (uint64_t*)calloc(Hz1.p, Hz1.p * sizeof(uint64_t));

	if (!msb_step1(msb_Hz1, msb_inf_Hz1, msb_sup_Hz1, &Hz1, wl))
	{
		fprintf(stderr, "Could not determine the MSB for step 1. \n");
		filter_free(&Hz1);
		free(msb_inf_Hz1);
		free(msb_sup_Hz1);
		freeMPFIMatrix(msb_Hz1,Hz1.p, 1);
		return -1;
	}	
	
	  //fprintf(stderr, "\n Step 1: MSB in mpfi format\n");
	  //writeMPFIMatrix(stderr, msb_Hz1, Hz1.p, 1, 5);
	  //printf("Step 1: MSB interval infinum \n");
	  //UINT64MatrixPrint(stderr, msb_inf_Hz1, 1, Hz1.p);
	  //printf("\n Step 1:  interval supremum \n");
	  //UINT64MatrixPrint(stderr, msb_sup_Hz1,  1, Hz1.p);

	 // fprintf(stderr, "------------------------------ Computing MSB on step 2 --------------------- \n");

	filter Delta;
	filter_allocate(&Delta, H->n, H->p + H->n, H->p + H->n);

	// first, we try to set MSB to inferior part of MSB interval from step 1
	createFilterDelta(&Delta, H, msb_inf_Hz1, wl);	

	mpfi_t *msb_step2_interval;
	msb_step2_interval = allocateMPFIMatrix(Hz1.p, 1, 64);
	uint64_t *msb_inf_step2;
	msb_inf_step2 = (uint64_t*)malloc(Hz1.p * sizeof(uint64_t)); 
	uint64_t *msb_sup_step2;
	msb_sup_step2 = (uint64_t*)malloc(Hz1.p * sizeof(uint64_t)); 

	/* We set the LSB constraint to the supremum of the MSB interval, obtained on the step 1.
		Then, the condition of breaking from the refinement loop is successful MSB check, 
		or if the LSB after refinment became larger than initially estimated MSB. 
	*/
	uint64_t *lsb_max;
	lsb_max = (uint64_t*)calloc(Hz1.p, Hz1.p * sizeof(uint64_t));
	UINT64MatrixCopy(lsb_max, msb_sup_Hz1, 1, Hz1.p);


	if (!msb_step2(msb_step2_interval, msb_inf_step2, msb_sup_step2, &Hz1, &Delta, wl,lsb_max, result))
	{
		fprintf(stderr, "Could not determine the MSB for step 2. \n");
		filter_free(&Hz1);
		free(msb_inf_Hz1);
		free(msb_sup_Hz1);
		freeMPFIMatrix(msb_Hz1,Hz1.p, 1);
		filter_free(&Delta);
		free(msb_inf_step2);
		free(msb_sup_step2);
		freeMPFIMatrix(msb_step2_interval,Hz1.p, 1);
		free(lsb_max);
		return -1;
	}

	
	//fprintf(stderr, "\n Step 2: MSB in mpfi format\n");
	//writeMPFIMatrix(stderr, msb_step2_interval, Hz1.p, 1, 5);
	//printf("Step 2: MSB \n");
	//UINT64MatrixPrint(stderr, msb_sup_step2, 1, Hz1.p);

	uint64_t *msb_new = (uint64_t*)calloc(Hz1.p, Hz1.p * sizeof(uint64_t)); 
	uint64_t *msb_z = (uint64_t*)calloc(Hz1.p, Hz1.p * sizeof(uint64_t)); 

	int flag_OK = 0;
	int steps = 0;
	if(!checkMSB(msb_new, msb_inf_Hz1, msb_sup_Hz1, msb_sup_step2, 1, Hz1.p))	//if format changed 
	{
		/*We are in the case when first check of MSB has not passed. */
		 //fprintf(stderr, "\n <-----Result----> The first check if msb_step2 is equal has not passed.\n");

		while(!flag_OK && steps < 10)
		{

			 //fprintf(stderr, "We are in the case, when msb_step2 is larger than msb_Hz1. \n");
			 //fprintf(stderr, "Steps = %d: MSB msb_step2_interval in mpfi format\n", steps );
			 //writeMPFIMatrix(stderr, msb_step2_interval, Hz1.p, 1, 25);
			 //fprintf(stderr, "Steps = %d: MSB msb_Hz1 in mpfi format from Step 1. \n", steps );
			 //writeMPFIMatrix(stderr, msb_Hz1, Hz1.p, 1, 25);

			steps++;
			/* Set the current estimation on MSBs to the corrected vector.
			 Update the LSBs in filter Delta after MSB vector correction. */
			UINT64MatrixCopy(msb_z, msb_new, 1, Hz1.p);
			filter_set_ur_to_eps(&Delta, msb_z, wl);	

			/* Try again Step 2. */
			if (!msb_step2(msb_step2_interval, msb_inf_step2, msb_sup_step2, &Hz1, &Delta, wl, lsb_max,result))
			{
				fprintf(stderr, "ERROR: Could not determine the MSB for step %d. \n", 2 + steps);
				free(msb_new);
				free(msb_z);
				filter_free(&Hz1);
				free(msb_inf_Hz1);
				free(msb_sup_Hz1);
				freeMPFIMatrix(msb_Hz1,Hz1.p, 1);
				filter_free(&Delta);
				free(msb_inf_step2);
				free(msb_sup_step2);
				freeMPFIMatrix(msb_step2_interval,Hz1.p, 1);
				free(lsb_max);
				return -1;
			}
				 // printf("Computing MSB on step %d -----------------\n",2 + steps  );
				 // fprintf(stderr, "Step %d: MSB:\n", 2 + steps );
				 // writeMPFIMatrix(stderr, msb_step2_interval, 1, Hz1.p, 5);
				 // printf("Step %d : MSB in mpfi format\n", 2 + steps  );
				 // UINT64MatrixPrint(stderr, msb_inf_step2, 1, Hz1.p);
				  // printf("Step %d : MSB in mpfi format\n", 2 + steps  );
				   //UINT64MatrixPrint(stderr, msb_sup_step2,  1, Hz1.p);
				   //fprintf(stderr, "Step %d: LSB:\n", 2 + steps );

				   	for(int i = 0; i < Hz1.p; ++i)
		            {
			            //fprintf(stderr," %d \t ", msb_sup_step2[i] - wl[i] + 1);
			            //fprintf(stderr, "compute LSB: lsb[%d] = msb[i] - w[i] + 1 = %lld - %lld + 1 = \n", i, msb_sup_step2[i], wl[i], msb_sup_step2[i] - wl[i] + 1 );
		            }
		            fprintf(stderr, "\n");

			/* Check if the new MSBs coincide with the current estimation. 
			If they do, we step out of the loop. Otherwise, correct the MSB vector. */
			if(checkMSB(msb_new, msb_z, msb_z, msb_sup_step2, 2, Hz1.p))
			{
				flag_OK = 1;
			}

		}
	}
	/* If the MSB positions computed from the first iteration of the 
	Step 2 passed the check, set the flag to 1. */
	else flag_OK = 1;

	end=clock();
		double overalltime = (double)(end - begin) / CLOCKS_PER_SEC;
		//fprintf(stderr, "\n\n ->Overall time: %f \n", overalltime);

	if(flag_OK)
	{
		int i;
		for( i = 0; i < Hz1.p; ++i)
		{
			lsb[i] = msb_sup_step2[i] - wl[i] + 1;
			msb[i] = msb_sup_step2[i];
		}

		//PerformWCPGSimulation(&Hz1, wl, msb );
		//PerformWCPGSimulation_forY(H, wl, msb );

		free(msb_new);
		free(msb_z);
		filter_free(&Hz1);
		free(msb_inf_Hz1);
		free(msb_sup_Hz1);
		freeMPFIMatrix(msb_Hz1,Hz1.p, 1);
		filter_free(&Delta);
		free(msb_inf_step2);
		free(msb_sup_step2);
		freeMPFIMatrix(msb_step2_interval,Hz1.p, 1);
		free(lsb_max);


		fxpf_result_setAdditionalSteps(result, steps);
		fxpf_result_setW(result, wl);
		fxpf_result_setLSB(result, lsb);
		fxpf_result_setMSB(result, msb);

		return steps;
	}
	else
	{
		fprintf(stderr, "ERROR: Could not determine the FxPF.\n");
		free(msb_new);
		free(msb_z);
		filter_free(&Hz1);
		free(msb_inf_Hz1);
		free(msb_sup_Hz1);
		freeMPFIMatrix(msb_Hz1,Hz1.p, 1);
		filter_free(&Delta);
		free(msb_inf_step2);
		free(msb_sup_step2);
		freeMPFIMatrix(msb_step2_interval,Hz1.p, 1);
		free(lsb_max);
		return -1;
	}
}

/* Checks the MSB vectors.
Returns 1 if OK, and zero if some changes were applied to vector msb_new */
int checkMSB(uint64_t *msb_new, uint64_t *msb1_inf, uint64_t *msb1_sup, uint64_t *msb2, uint64_t stage, uint64_t length)
{
	// fprintf(stderr, "\nChecking MSB_diamond against MSB....\n");
	// fprintf(stderr, "\nMSB_z:\n");
	// UINT64MatrixPrint(stderr, msb1_sup, 1, length);
	// fprintf(stderr, "\nMSB_z_diamond:\n");
	 // UINT64MatrixPrint(stderr, msb2, 1, length);
	
	int i;
	int flag = 1;
	for(i = 0; i < length; ++i)
	{
		if(stage == 1)
		{
			// printf("msb2[%d] = %llu, msb1_inf[%d] = %llu\n", i,msb2[i], i, msb1_inf[i] );
			if(msb2[i] > msb1_inf[i])	//MSB computed on step2 is larger than  msb1_inf
			{
				flag = 0;
				msb_new[i] = msb2[i];
			}
			else
				msb_new[i] = msb1_inf[i];
		}
		else if (stage == 2)
		{
		 //printf("msb2[%d] = %llu, msb1_sup[%d] = %llu\n", i,msb2[i], i, msb1_sup[i] );
			if(msb2[i] > msb1_sup[i])
			{
				flag = 0;
				msb_new[i] = msb1_sup[i] + 1;
			}
			else
				msb_new[i] = msb2[i];
		}
	}

	return flag;
}

void PerformWCPGSimulation_forY(filter *Hz1, uint64_t *wl, uint64_t *msb )
{

	/*-------------------- */
			mpfr_prec_t prec = 64;

			uint64_t n = Hz1->n;
			uint64_t p = Hz1->p;
			uint64_t q = Hz1->q;

			mpfr_t *A;
			A = allocateMPFRMatrix(n, n, prec);
			doubleToMPFRMatrix(A, Hz1->A, n, n);

			mpfr_t *B;
			B = allocateMPFRMatrix(n, q, prec);
			doubleToMPFRMatrix(B, Hz1->B, n, q);

			mpfr_t *C;
			C = allocateMPFRMatrix(p, n, prec);
			doubleToMPFRMatrix(C, Hz1->C, p, n);

			mpfr_t *D;
			D = allocateMPFRMatrix(p, q, prec);
			doubleToMPFRMatrix(D, Hz1->D, p, q);


			int *lsb_x = (int*)malloc(n * sizeof(int));
			int *lsb_y = (int*)malloc(p * sizeof(int));
			int *msb_x = (int*)malloc(n * sizeof(int));
			int *msb_y = (int*)malloc(p * sizeof(int));
			int i;
			for(i = 0; i < n; ++i)
			{
				lsb_x[i] = msb[i] - wl[i] + 1;
				msb_x[i] = msb[i];
				//printf("msb_x[i] = %d\n", msb_x[i]);
			}
			for(i = 0; i < p; ++i)
			{
				lsb_y[i] = msb[n + i] - wl[n + i] + 1;
				msb_y[i] = msb[n + i];
				//printf("msb_y[i] = %d\n", msb_y[i]);
			}


			FILE *input = fopen("./examples/RAIM_inputWCPG.txt", "r");
			FILE *output_exact = fopen("./examples/RAIM1_output_exact.txt", "w");
			FILE *output_q = fopen("./examples/RAIM1_output_quantized.txt", "w");

			if(!input)
				fprintf(stderr, "Error opening file. \n");
			if(!output_exact)
				fprintf(stderr, "Error opening file. \n");
			if(!output_q)
				fprintf(stderr, "Error opening file. \n");

			int Tfinal;
			int tt;

			mpfr_t *U_i;
			fscanf(input, "%d %d<", &Tfinal, &tt);
			//Tfinal = 11;
			U_i = allocateMPFRMatrix(Tfinal, Hz1->q, 64);
			printf("====> Simulating the output with Tfinal = %d \n", Tfinal);
			//printf("simulating the output %d, and %d \n", -n, p);


			mpfr_t *X_out = allocateMPFRMatrix(n, Tfinal, prec);
			mpfr_t *Y_out = allocateMPFRMatrix(p, Tfinal, prec);
	/*
			printf("Simulating the system\n");
			printf("Matrix A\n");
			writeMPFRMatrix(stderr, A, n, n, 5, MPFR_RNDN);
			printf("Matrix B\n");
			writeMPFRMatrix(stderr, B, n, q, 5, MPFR_RNDN);
			printf("Matrix C\n");
			writeMPFRMatrix(stderr, C, p, n, 5, MPFR_RNDN);
			printf("Matrix D\n");
			writeMPFRMatrix(stderr, D, p, q, 5, MPFR_RNDN);
	*/
			fprintf(output_exact, "%d %d %d", n, p, Tfinal);
			fprintf(output_q, "%d %d %d", n, p, Tfinal);

			for(i = n+p-1; i < n + p; ++i)
			{
				 int index = i - n;

				 //int index = 0;
				 readMPFRMatrix(U_i, input, Tfinal, tt, MPFR_RNDN);
				// writeMPFRMatrix(stderr, U_i, tt, Tfinal, 10, MPFR_RNDN);
				 //printf("simulating the output \n");
				 //fprintf(output_exact, "Output exact \n");
				 //fprintf(output_q, "Output quantized  \n");
				 fprintf(output_exact, "\n%d 1\n", Tfinal);
				 fprintf(output_q, "\n%d 1\n", Tfinal);
				 simulateSS_SIMO(output_exact, output_q, index, X_out, Y_out, A, B, C, D, U_i, n, p, q, Tfinal, lsb_x, lsb_y, msb_x, msb_y);
				 //fprintf(output_exact, "];\n");
				 //fprintf(output_q, "];\n");
				 //printf("Verdict of the simulation %d: %d \n", i, simResult);
				 index++;
			}

}



void PerformWCPGSimulation(filter *Hz1, uint64_t *wl, uint64_t *msb )
{

	/*-------------------- */
			mpfr_prec_t prec = 64;

			uint64_t n = Hz1->n;
			uint64_t p = Hz1->p;
			uint64_t q = Hz1->q;

			mpfr_t *A;
			A = allocateMPFRMatrix(n, n, prec);
			doubleToMPFRMatrix(A, Hz1->A, n, n);

			mpfr_t *B;
			B = allocateMPFRMatrix(n, q, prec);
			doubleToMPFRMatrix(B, Hz1->B, n, q);

			mpfr_t *C;
			C = allocateMPFRMatrix(p, n, prec);
			doubleToMPFRMatrix(C, Hz1->C, p, n);

			mpfr_t *D;
			D = allocateMPFRMatrix(p, q, prec);
			doubleToMPFRMatrix(D, Hz1->D, p, q);


			int *lsb_x = (int*)malloc(n * sizeof(int));
			int *lsb_y = (int*)malloc(p * sizeof(int));
			int *msb_x = (int*)malloc(n * sizeof(int));
			int *msb_y = (int*)malloc(p * sizeof(int));
			int i;
			for(i = 0; i < n; ++i)
			{
				lsb_x[i] = msb[i] - wl[i] + 1;
				msb_x[i] = msb[i];
			}
			for(i = 0; i < p; ++i)
			{
				lsb_y[i] = msb[i] - wl[i] + 1;
				msb_y[i] = msb[i];
			}




			FILE *input = fopen("./examples/RAIM_inputWCPG.txt", "r");
			FILE *output_exact = fopen("./examples/RAIM1_output_exact.txt", "w");
			FILE *output_q = fopen("./examples/RAIM1_output_quantized.txt", "w");

			if(!input)
				fprintf(stderr, "Error opening file. \n");
			if(!output_exact)
				fprintf(stderr, "Error opening file. \n");
			if(!output_q)
				fprintf(stderr, "Error opening file. \n");

			int Tfinal;
			int tt;

			mpfr_t *U_i;
			fscanf(input, "%d %d", &Tfinal, &tt);
			//Tfinal = 11;
			U_i = allocateMPFRMatrix(Tfinal, Hz1->q, 64);
			printf("====> Simulating the output with Tfinal = %d \n", Tfinal);
			//printf("simulating the output %d, and %d \n", -n, p);


			mpfr_t *X_out = allocateMPFRMatrix(n, Tfinal, prec);
			mpfr_t *Y_out = allocateMPFRMatrix(p, Tfinal, prec);
	/*
			printf("Simulating the system\n");
			printf("Matrix A\n");
			writeMPFRMatrix(stderr, A, n, n, 5, MPFR_RNDN);
			printf("Matrix B\n");
			writeMPFRMatrix(stderr, B, n, q, 5, MPFR_RNDN);
			printf("Matrix C\n");
			writeMPFRMatrix(stderr, C, p, n, 5, MPFR_RNDN);
			printf("Matrix D\n");
			writeMPFRMatrix(stderr, D, p, q, 5, MPFR_RNDN);
	*/
			fprintf(output_exact, "%d %d %d", n, p, Tfinal);
			fprintf(output_q, "%d %d %d", n, p, Tfinal);

			for(i = n+p-1; i < n + p; ++i)
			{
				 int index = i - n;

				 //int index = 0;
				 readMPFRMatrix(U_i, input, Tfinal, tt, MPFR_RNDN);
				// writeMPFRMatrix(stderr, U_i, tt, Tfinal, 10, MPFR_RNDN);
				 //printf("simulating the output \n");
				 //fprintf(output_exact, "Output exact \n");
				 //fprintf(output_q, "Output quantized  \n");
				 fprintf(output_exact, "\n%d 1\n", Tfinal);
				 fprintf(output_q, "\n%d 1\n", Tfinal);
				 simulateSS_SIMO(output_exact, output_q, index, X_out, Y_out, A, B, C, D, U_i, n, p, q, Tfinal, lsb_x, lsb_y, msb_x, msb_y);
				 //fprintf(output_exact, "];\n");
				 //fprintf(output_q, "];\n");
				 //printf("Verdict of the simulation %d: %d \n", i, simResult);
				 index++;
			}

}






/*
	For a LTI filter in State-Space representation H=(A,B,C,D),
	with n states, q inputs, p outputs;
	and an input vector U(k) of size q x T in form 
		Y(k) = C * x(k) + D * u(k)
		X(k+1) = A * x(k) + B * u(k)
		X(0) = 0 
	the function simulates the implementation of the filter
	in Fixed-Point Format given in the argument.

	At each time instance the function computes the state and output exactly 
	and then performs correct rounding to the Fixed-Point format.

	The function returns a boolean indicating wether an overflow occured during the simulation.

	Input:
		mpfr_t *X_out - T x n MPFR matrix, wich on the output will contain the quantized states for k \in [0;T]
		mpfr_t *Y_out - T x p MPFR matrix, wich on the output will contain the quantized outputs for k \in [0;T]
		mpfr_t *A 	  - n x n MPFR matrix defining the filter
		mpfr_t *B 	  - n x q MPFR matrix defining the filter
		mpfr_t *C 	  - p x n MPFR matrix defining the filter
		mpfr_t *D 	  - p x q MPFR matrix defining the filter
		mpfr_t *U 	  - T x q MPFR matrix, each line U[k,:] is a 1xq vector containing the input at time k
		uint64_t n 	  - number of states
		uint64_t p 	  - number of outputs
		uint64_t q 	  - number of inputs
		uint64_t T 	  - simulation length
		int *lsb_x 	  - 1 x n vector containing the LSB positions for state vector implementation
		int *lsb_y 	  - 1 x q vector containing the LSB positions for output vector implementation
		int *msb_x 	  - 1 x n vector containing the MSB positions for state vector implementation
		int *msb_y    - 1 x p vector containing the MSB positions for output vector implementation
 */


int simulateSS_SIMO(FILE *out_exact, FILE *out_q, int index, mpfr_t *X_out, mpfr_t *Y_out, mpfr_t *A, mpfr_t *B, mpfr_t *C, mpfr_t *D, mpfr_t *U, uint64_t n, uint64_t p, uint64_t q, uint64_t T, int *lsb_x, int*lsb_y, int *msb_x, int *msb_y)
{


	//printf("SImulating output. \n");

	int flag = 1;
	mpfr_prec_t prec = 64;

	mpfr_t *X;
	X = allocateMPFRMatrix(1, n, prec);
	setMatrixZero(X, n, 1);

	mpfr_t *Y;
	Y= allocateMPFRMatrix(1, p, prec);
	setMatrixZero(Y, p, 1);

	mpfr_t *X_exact;
	X_exact = allocateMPFRMatrix(1, n, prec);
	setMatrixZero(X_exact, n, 1);

	mpfr_t *Y_exact;
	Y_exact = allocateMPFRMatrix(1, p, prec);
	setMatrixZero(Y_exact, p, 1);

	mpfr_t *X_q;
	X_q = allocateMPFRMatrix(1, n, prec);
	setMatrixZero(X_q, n, 1);

	mpfr_t *Y_q;
	Y_q = allocateMPFRMatrix(1, p, prec);
	setMatrixZero(Y_q, p, 1);

	int i, k;

	mpfr_t *u_k;
	u_k = allocateMPFRMatrix(q, 1, prec);

	mpfr_t *t1, *t2;
	t1 = allocateMPFRMatrix(p, 1, prec);
	t2 = allocateMPFRMatrix(p, 1, prec);
	mpfr_t *t3, *t4;
	t3 = allocateMPFRMatrix(n, 1, prec);
	t4 = allocateMPFRMatrix(n, 1, prec);

	mpfr_t bound;
	mpfr_init(bound);

	mpfr_t tmp;
	mpfr_init(tmp);

	/*
	if(index >= 0)
	{
		fprintf(out_exact, "Y%d = [ \n", index);
		fprintf(out_q, "Y_quantized%d = [\n", index);
	}
	else
	{
		fprintf(out_exact, "X%llu = [ \n", index + n);
		fprintf(out_q, "X_quantized%llu = [\n", index + n);
	}
*/

	mpfr_t scratch;
	mpfr_init(scratch);

	mpfr_t *scratch2 = allocateMPFRMatrix(2,1,64);

	for(k = 0; k < T; ++k)
	{
		printf("%d\n", k);




		/* We first extract the input vector for the time k. */	
		MatrixExtractRow(u_k, U, T, q, k);		
		
		/* Compute Y_exact as exact value of SoP C * x_exact + D * u(k, :) 
		We compute:
		u_k = u[k,:];
		t1 = D * u_k
		t2 = C * x_exact 
		Y_exact = t1 + t2
		*/
		MatrixMulVector_exact(t1, D, u_k, p, q, scratch2);
		MatrixMulVector_exact(t2, C, X_exact, p, n, scratch2);
		MatrixAdd_exact(Y_exact, t1, t2, p, 1);


	
		/* We quantize Y_exact */

		VectorQuantize(Y_q, Y_exact, lsb_y, p, scratch);
	
		/* Now we compute the next X_exact:
			X_exact = A * x_exact + B * u(k,:)
		We compute:
		u_k = u[k,:];
		t3 = B * u_k
		t4 = A * x_exact
		X_exact = t1 + t2 */
		MatrixMulVector_exact(t3, B, u_k, n, q, scratch2);
		MatrixMulVector_exact(t4, A, X_exact, n, n, scratch2);
		MatrixAdd_exact(X_exact, t3, t4, n, 1);
	
		/* We quantize the exact state vector */
		VectorQuantize(X_q, X_exact, lsb_x, n, scratch);





		/* Set the X_exact to the quantized value */
		MPFRMatrixCopy(X_exact, X_q, 1, n);

		/* Saving the values of the quantized output to the vector X_out, Y_out */
		MatrixSetRow(X_out, X_q, T, n, k);
		MatrixSetRow(Y_out, Y_q, T, p, k);


		//Writing the output with index 'index' to the file output_q
		// mpfr_out_str(out_q, 10, 10, Y_q[index], MPFR_RNDN);
		// fprintf(out_q, "\n");


		if(index >= 0)
		{
			mpfr_out_str(out_q, 10, 20, Y_q[index], MPFR_RNDN);
			fprintf(out_q, "\n");
		}
		else
		{
			mpfr_out_str(out_q, 10, 20, X_q[n + index], MPFR_RNDN);
			fprintf(out_q, "\n");

		}
	

		/*------------------------------------------------
			COMPUTING THE EXACT OUTPUT AND STATE 
			and writing them to the out_exact file
		--------------------------------------------------*/

		MatrixMulVector_exact(t1, D, u_k, p, q, scratch2);
		MatrixMulVector_exact(t2, C, X, p, n, scratch2);
		MatrixAdd_exact(Y, t1, t2, p, 1);

		//printf("===> y_exact for k=%d\n", k);
		//writeMPFRMatrix(stderr, Y, 1, p, 10, MPFR_RNDN);
	//	printf("===> y_quantized for k=%d\n", k);
	//	writeMPFRMatrix(stderr, Y_q, 1, p, 10, MPFR_RNDN);

		MatrixMulVector_exact(t3, B, u_k, n, q, scratch2);
		MatrixMulVector_exact(t4, A, X, n, n, scratch2);
		MatrixAdd_exact(X, t3, t4, n, 1);

		//Writing the output with index 'index' to the file output_EXACT
		if(index >= 0)
		{
			
			mpfr_out_str(out_exact, 10, 20, Y[index], MPFR_RNDN);
			fprintf(out_exact, "\n");
		}
		else
		{
			mpfr_out_str(out_exact, 10, 20, X[n + index], MPFR_RNDN);
			fprintf(out_exact, "\n");

		}
		


		/* We check each quantized output. If 
			abs(Y_q[i]) >= 2^msb 
			then we have an overflow.
			Note: actually if Y_q[i] == -2^msb then there is no overflow
			but we, in our initial conditions have asked |Y_q| <= 2^msb - 2^lsb.

			 */
		
		for(i = 0; i < p; ++i)
		{
			
			mpfr_exp_t msb_i = msb_y[i];
			mpfr_set_ui_2exp(bound, 1, msb_i, MPFR_RNDN);		//bound = 2^msb

			mp_prec_t precY = mpfr_get_prec(Y_q[i]);
			mpfr_set_prec(tmp, precY);
			mpfr_abs(tmp, Y_q[i], MPFR_RNDN);		//tmp = abs(Y_q[i])

			if(mpfr_cmp(tmp, bound )>= 0)
			{
				fprintf(stderr, "Warning: overflow during simulation ! Exists a value of quantized output Y_q such that it is larger than bound.\n");
				fprintf(stderr,"The bound is: \n");
				mpfr_out_str(stderr, 10, 10, bound, MPFR_RNDD);
				fprintf(stderr,"\n The |Y[%d]| at k=%d value is: \n", i, k);
				mpfr_out_str(stderr, 10, 10, tmp, MPFR_RNDN);
				flag = 0;
			}

		}

		/* We check each quantized output. If 
			abs(X_q[i]) >= 2^msb 
			then we have an overflow.
			Note: actually if X_q[i] == -2^msb then there is no overflow
			but we, in our initial conditions have asked |X_q| <= 2^msb - 2^lsb.

			 */
			
		for(i = 0; i < n; ++i)
		{
			
			mpfr_exp_t msb_i = msb_x[i];
			mpfr_set_ui_2exp(bound, 1, msb_i, MPFR_RNDN);		//bound = 2^msb - 2^lsb

			mp_prec_t precX = mpfr_get_prec(X_q[i]);
			mpfr_set_prec(tmp, precX);
			mpfr_abs(tmp, X_q[i], MPFR_RNDN);					//tmp = abs(X_q[i])

			if(mpfr_cmp(tmp, bound )>= 0)
			{
				fprintf(stderr, "Warning: overflow during simulation ! Exists a value of quantized state X_q such that it is larger than bound.\n");
				printf("\n The bound is: \n");
				mpfr_out_str(stderr, 10, 10, bound, MPFR_RNDD);
				printf("\n The |X[%d]| at k=%d value is: \n", i, k);
				mpfr_out_str(stderr, 10, 10, tmp, MPFR_RNDN);
				flag = 0;
			}

			
		}
			

		
	}
		


		mpfr_clear(bound);
		mpfr_clear(tmp);
		freeMPFRMatrix(t1, 1, p);
		freeMPFRMatrix(t2, 1, p);
		freeMPFRMatrix(t3, 1, n);
		freeMPFRMatrix(t4, 1, n);
		freeMPFRMatrix(u_k, 1, q);
		freeMPFRMatrix(X_q, 1, n);
		freeMPFRMatrix(X_exact, 1, n);
		freeMPFRMatrix(Y_q, 1, p);
		freeMPFRMatrix(Y_exact, 1, p);
		freeMPFRMatrix(Y, 1, p);
		freeMPFRMatrix(X, 1, n);
	
		return flag;


}


/* A wrapper function that computes the Fixed-Point Formats for a given filter that is defined via State-Space matrices.
 * This function calls determineFXPF with appropriate parameters and then parses the results.
 * This function is present for easier data transfer to Python using ctypes.
 *
 * All output arguments (msb, lsb, error, additionalSteps)
 */
int FXPF(uint64_t *msb, int *lsb, int *additionalSteps, double *error,
int *wl, double *A, double *B, double *C, double *D, int n, int p, int q, double *u_bound)
{

    filter H;

    filter_allocate(&H, n, p, q);
	filter_set(&H, n, p, q, A, B, C, D, u_bound);

	fxpf_result result;
	fxpf_result_allocate(&result, H.p + H.n);


    if(determineFXPF(msb, lsb, &H, wl, &result) < 0)
	{
		fprintf(stderr, "Error determining Fixed-Point Formats.\n");
		return -1;
	}
    else
    {

        *additionalSteps = result.additionalSteps;

	    for(int i = 0; i < n+p; ++i)
	    {
		    error[i] = mpfi_get_d(result.error[i]);
	    }
/*
	        fprintf(stderr, "####################################################################\n");
	        fprintf(stderr, "The filter:\n");
		    filter_print(stderr, &H);

            fprintf(stderr, "\n\n\n--------- FXPF Result ----------\n");
		    fprintf(stderr, "Additional steps: %d,\n", *additionalSteps);
		    fprintf(stderr, "--------------------------------------------------------------\n");
		    fprintf(stderr, "Wordlengths: \n");
		    UINT64MatrixPrint(stderr, wl, 1, p + n);
		    fprintf(stderr, "MSBs: \n");
		    UINT64MatrixPrint(stderr, msb, 1, p + n);
		    fprintf(stderr, "LSBs: \n");
		    INTMatrixPrint(stderr, lsb, 1, p + n);
		    fprintf(stderr, "Errors: \n");
		    DoubleMatrixPrint(stderr, error, 1, p + n);
		    fprintf(stderr, "####################################################################\n");
*/
	}

    filter_free(&H);
	fxpf_result_free(&result);

	return 0;
}















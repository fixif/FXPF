/*
 ============================================================================
 Name        : FXPF.h
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : 	  

 ============================================================================
*/

#ifndef _FXPF_H_
#define _FXPF_H_


#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <mpfi.h>
#include <mpfi_io.h>
#include <wcpg.h>
#include <limits.h>

#include "matrixalg.h"
#include "mpfi_matrixalg.h"
#include "filter.h"


typedef struct fxpf_result_struct
{
	int additionalSteps;
	uint64_t length;
	mpfi_t *error;
	uint64_t *msb;
	int *lsb;
	uint64_t *w;
}fxpf_result;




void fxpf_result_allocate					(fxpf_result *result, uint64_t length);
void fxpf_result_free						(fxpf_result *result);
void fxpf_result_setError					(fxpf_result *result, mpfi_t *error);
void fxpf_result_setAdditionalSteps			(fxpf_result *result, int steps);
void fxpf_result_setW						(fxpf_result *result, uint64_t *w);
void fxpf_result_setMSB						(fxpf_result *result, uint64_t *msb);
void fxpf_result_setLSB						(fxpf_result *result, int *lsb);
void fxpf_result_print						(FILE *file, fxpf_result *result);






/* Computes the precision of WCPG computation of filter H such that
the MSB computed with above mentioned WCPG is either exact, either overestimated by one.
The precision eps is computed via following formula:
	eps = max_i {eps_i}
	eps_i <= (WCPGLowerBound(H) * ur)_i / alpha * sum_j {ur_j} for i=1:p

with alpha = 1 or 2.
Note: we are going to compute WCPG lower bound with precision 64.
 */
int getWCPGeps								(mpfr_t epsWCPG, filter *H, uint64_t alpha);

int msb_step1								(mpfi_t *msb, uint64_t *msb_left, uint64_t *msb_right,
											 filter *H, uint64_t *w);

int msb_step2								(mpfi_t *msb, uint64_t *msb_left, uint64_t *msb_right,
											 filter *Hz1, filter *Delta, uint64_t *w, uint64_t *lsb_max,
											 fxpf_result *result);

int determineFXPF							(uint64_t *msb, int *lsb, filter *H, uint64_t *wl, fxpf_result *result);

int FXPF                                    (uint64_t *msb, int *lsb, int *additionalSteps, double *error, int *wl,
                                            double *A, double *B, double *C, double *D,
                                            int n, int p, int q, double *u_bound);



int computeLSB								(int *lsb, uint64_t *msb, uint64_t *w, uint64_t length);
int checkMSB								(uint64_t *msb_new, uint64_t *msb1_inf, uint64_t *msb1_sup,
											 uint64_t *msb2, uint64_t stage, uint64_t length);

void PerformWCPGSimulation_forY				(filter *Hz1, uint64_t *wl, uint64_t *msb );
void PerformWCPGSimulation					(filter *Hz1, uint64_t *wl, uint64_t *msb );
int simulateSS_SIMO							(FILE *out_exact, FILE *out_q, int index,
											 mpfr_t *X_out, mpfr_t *Y_out, mpfr_t *A,
											 mpfr_t *B, mpfr_t *C, mpfr_t *D, mpfr_t *U,
											 uint64_t n, uint64_t p, uint64_t q,
											 uint64_t T, int *lsb_x, int*lsb_y,
											 int *msb_x, int *msb_y);

#endif

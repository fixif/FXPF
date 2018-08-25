/*
 ============================================================================
 Name        : filter.c
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : 	  

 ============================================================================
*/

#include "filter.h"

void filter_allocate(filter *f, uint64_t n, uint64_t p, uint64_t q)
{
	f->n = n;
	f->p = p;
	f->q = q;
	f->A = (double*)malloc(n * n * sizeof(double));
	f->B = (double*)malloc(n * q * sizeof(double));
	f->C = (double*)malloc(p * n * sizeof(double));
	f->D = (double*)malloc(p * q * sizeof(double));
	f->ur = (double*)malloc(1 * q * sizeof(double));
	f->w = (uint64_t*)malloc(1 * p * sizeof(uint64_t));
	f->msbU = (int*)malloc(1 * q * sizeof(int));
	f->lsbU = (int*)malloc(1 * q * sizeof(int));
}

void filter_free(filter *f)
{
	free(f->A);
	free(f->B);
	free(f->C);
	free(f->D);
	free(f->ur);
	free(f->w);
	free(f->msbU);
	free(f->lsbU);
}

void filter_set(filter *f, uint64_t n, uint64_t p, uint64_t q, double *A, \
					 double *B, double *C, double *D, double *ur)
{
	if(n != f->n || p != f->p || q != f->q )
		fprintf(stderr, "ERROR: Trying to initilize filter struct with wrong sizes \n");


	DoubleMatrixCopy(f->A, A, n, n);
	DoubleMatrixCopy(f->B, B, n, q);
	DoubleMatrixCopy(f->C, C, p, n);
	DoubleMatrixCopy(f->D, D, p, q);
	DoubleMatrixCopy(f->ur, ur, 1, q);

}

void filter_set_uFormat(filter *f, int *msb, int *lsb)
{
	int i;
	for(i = 0; i < f->q; ++i)
	{
		f->msbU[i] = msb[i];
		f->lsbU[i] = lsb[i];
	}
}

void filter_set_w(filter *f, uint64_t *w_z)
{
	int i;
	for(i = 0; i < f->p; ++i)
	{
		f->w[i] = w_z[i];
	}
}
void filter_set_ur_to_eps(filter *f, uint64_t *msb, uint64_t *w)
{
	int tmp;
	int i;
	for(i = 0; i < f->q; ++i)
	{
			tmp = (int)(msb[i] - w[i] + 1);
			//printf("msb[%d] = %llu, w[%d] = %llu \n", i, msb[i], i, w[i]);
			f->ur[i] = pow(2, tmp);
			//newU[i] = one >> tmp;
			//printf("newU = %f, newU = %e \n",newU[i], newU[i]);
	}
}
void filter_set_w_eql(uint64_t* wl, uint64_t w_z, uint64_t length)
{
	int i;
	for(i = 0; i < length; ++i)
	{
		wl[i] = w_z;
	}
}

void filter_print(FILE *file, filter *f)
{
	fprintf(file, "Filter order: %llu.\nNumber of inputs: %llu \nNumber of outputs: %llu \n",f->n, f->q, f->p);
	fprintf(file, "Input radius ur: \n");
	DoubleMatrixPrint(file, f->ur, 1, f->q);
	fprintf(file, "Matrix A: \n");
	DoubleMatrixPrint(file, f->A, f->n, f->n);
	fprintf(file, "Matrix B: \n");
	DoubleMatrixPrint(file, f->B, f->n, f->q);
	fprintf(file, "Matrix C: \n");
	DoubleMatrixPrint(file, f->C, f->p, f->n);
	fprintf(file, "Matrix D: \n");
	DoubleMatrixPrint(file, f->D, f->p, f->q);
}

/*
	Creates a filter Hz1, which is a following modification of filter H:
                  (I) (0)
	Hz1 := {A, B, (C),(D)}
	with output vector zeta = (x'(k) y'(k))' size (n + p) x 1.

	Space for Hz1 is assumed to be preallocated outside the function.
 */
void createFilterHz1(filter *Hz1, filter *H)
{
	double *newC, *newD;
	newC = (double*)malloc((H->n+H->p) * H->n * sizeof(double));
	newD = (double*)malloc((H->n+H->p) * H->q * sizeof(double));

	double *I, *Z;
	I = (double*)malloc(H->n * H->n * sizeof(double));
	Z = (double*)malloc(H->q * H->n * sizeof(double));
	DoubleMatrixIdent(I, H->n);
	DoubleMatrixZero(Z, H->n, H->q);

	doubleMatrixConcatVertical(newC, I, H->C, H->n, H->p, H->n);
	doubleMatrixConcatVertical(newD, Z, H->D, H->n, H->p, H->q);

	filter_set(Hz1, H->n, H->n + H->p, H->q, H->A, H->B, newC, newD, H->ur);

	free(I);
	free(Z);
	free(newC);
	free(newD);
}


/* Creates a filter Delta, which is the difference between impleneted and exact filters.
   This filter has following form:
                       (I) (0 0) 
	Delta := {A, (I 0),(C),(0 I)}

	with output vector of size (n + p) x 1                         (2^(LSB_x))
		input vector of size (n + p) x 1 for which the bound is ur=(2^(LSB_y)).

	LSB vectors are obtained from vector msb given as input argument
	via formula LSB = MSB - wordlength + 1.

	Space for Delta is assumed to be preallocated outside the function.

	Inputs:
		* filter H - the target for implementation, all fields specified
		* msb - vector of MSB for state- and output-variables of filter H,
				should be of size (H->n + H->p) x 1. 
		* wl - vector of wordlengths for ste- and output variables of filter H;
				should be of size (H->n + H->p) x 1. 
	Output:
		*filter Delta

 */

void createFilterDelta(filter *Delta, filter *H, uint64_t *msb, uint64_t *w)
{


	double *newB, *newC, *newD, *newU;
	newB = (double*)calloc(H->n * (H->n+H->p), H->n* (H->n+H->p)  * sizeof(double));
	newC = (double*)calloc((H->n+H->p) * H->n, (H->n+H->p) * H->n * sizeof(double));
	newD = (double*)calloc((H->n+H->p) * (H->n+H->p), (H->n+H->p) * (H->n+H->p) * sizeof(double));
	newU = (double*)calloc((H->n+H->p), (H->n+H->p) * sizeof(double));

	double *Inn, *Ipp, *Znp, *Znp_n, *ZI;
	Inn = (double*)calloc(H->n * H->n, H->n * H->n * sizeof(double));
	Ipp = (double*)calloc(H->p * H->p, H->p * H->p * sizeof(double));
	Znp = (double*)calloc(H->n * H->p, H->n * H->p * sizeof(double));
	Znp_n = (double*)calloc((H->n + H->p) * H->n, (H->n + H->p) * H->n * sizeof(double));
	ZI = (double*)calloc((H->n + H->p) * H->p , (H->n + H->p) * H->p * sizeof(double));
	DoubleMatrixIdent(Inn, H->n);
	DoubleMatrixIdent(Ipp, H->p);

	/* new matrix B */
	doubleMatrixConcatHorizontal(newB, Inn, Znp, H->n, H->n, H->p);	

	/* new Matrix C */
	doubleMatrixConcatVertical(newC, Inn, H->C, H->n, H->p, H->n);

	/* new Matrix D */
	doubleMatrixConcatVertical(ZI, Znp, Ipp, H->n, H->p, H->p);
	doubleMatrixConcatHorizontal(newD, Znp_n, ZI, (H->n + H->p), H->n, H->p);

	/* new Matrix ur */
	int tmp;
	int i;
	for(i = 0; i < H->q; ++i)
	{
			tmp = (int)(msb[i] - w[i] + 1);
			//printf("msb[%d] = %llu, w[%d] = %llu \n", i, msb[i], i, w[i]);
			newU[i] = pow(2, tmp);

			//newU[i] = one >> tmp;
			//printf("newU = %f, newU = %e \n",newU[i], newU[i]);
	}


	// fprintf(stderr, "Input radius for Delta filter: \n");
	// DoubleMatrixPrint(stderr, newU, 1, Delta->q);
	// fprintf(stderr, "Wordlengths for Delta filter: \n");
	// UINT64MatrixPrint(stderr, w, 1, Delta->p);

	filter_set(Delta, H->n, H->n + H->p, H->n + H->p, H->A, newB, newC, newD, newU);

	free(newB);
	free(newC);
	free(newD);
	free(newU);
	free(Inn);
	free(Ipp);
	free(Znp);
	free(Znp_n);
	free(ZI);

}


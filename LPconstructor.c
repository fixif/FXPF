/*
 ============================================================================
 Name        : LPconstructor.c
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : This file contains source files for the module responsible 
 				for the consutruction of the Integer Linear Programming problem
 				and writing it into the *.lp format. (which is later to be used 
 				with SCIP)

 ============================================================================
*/

#include "LPconstructor.h"


/* If x is NaN or Inf, return 0

   If x is 0, return 0

   Otherwise return an integer E such that

   2^-E * x is an odd integer.

*/
mp_exp_t exponentMpfr(mpfr_t x) {
  mp_exp_t E;
  mpz_t m;
  
  if (!mpfr_number_p(x)) return ((mp_exp_t) 0);
  if (mpfr_zero_p(x)) return ((mp_exp_t) 0);
  
  mpz_init(m);
  E = mpfr_get_z_2exp(m, x);
  mpz_abs(m,m);
  E += mpz_scan1(m, 0);
  mpz_clear(m);

  return E;
}


void minOnLine(mpfr_exp_t *v, mpfr_t *A, uint64_t n, uint64_t m)
{
	int i,j;
	mpfr_exp_t min_line;
	mpfr_exp_t current;
	for(i = 0; i < n; ++i)
	{
		min_line = exponentMpfr(A[i*m]);
		for(j = 0; j < m; ++j)
		{
			current = exponentMpfr(A[i*m + j]);
			if( current < min_line && current!=0 ) 				//<------Christoph's function
				min_line = exponentMpfr(A[i*m + j]);

		}
		v[i] = min_line;
	}
}


void LPcomputeRightSide(mpfr_t *rightSide, mpfr_t *A, mpfr_t *B, mpfr_t *C,mpfr_t *D,\
			 uint64_t n, uint64_t p, uint64_t q,\
			 mpfr_t *x_low, mpfr_t *x_up,\
			 mpfr_t *y_low, mpfr_t *y_up,\
			 mpfr_t *u_low, mpfr_t *u_up)
{
	/* We compute the right side vector */

	int i,j;
	// void MPFRMatrixConcatVertical(mpfr_t *C, mpfr_t *A, mpfr_t *B, uint64_t m, uint64_t k, uint64_t n);

	mpfr_t *Ax, *Bu, *Cx, *Du;
	Ax = allocateMPFRMatrix(n ,1, 64);
	Bu = allocateMPFRMatrix(n ,1, 64);
	Cx = allocateMPFRMatrix(p ,1, 64);
	Du = allocateMPFRMatrix(p ,1, 64);

	mpfr_t *scratch = allocateMPFRMatrix(2,1,64);

	MatrixMulVector_exact(Ax, A, x_low, n, n, scratch);
	MatrixMulVector_exact(Bu, B, u_low, n, q, scratch);
	MatrixMulVector_exact(Cx, C, x_low, p, n, scratch);
	MatrixMulVector_exact(Du, D, u_low, p, q, scratch);

	freeMPFRMatrix(scratch, 2,1);

	if(matrixIsNan_mpfr(Ax, n, 1))
	{
		fprintf(stderr, "Error: multiplying matrix A by state vector results in NaN matrix. \n");
		return;
	}
	if(matrixIsNan_mpfr(Bu, n, 1))
	{
		fprintf(stderr, "Error: multiplying matrix A by state vector results in NaN matrix. \n");
		return;
	}
	if(matrixIsNan_mpfr(Cx, p, 1))
	{
		fprintf(stderr, "Error: multiplying matrix A by state vector results in NaN matrix. \n");
		return;
	}
	if(matrixIsNan_mpfr(Du, p, 1))
	{
		fprintf(stderr, "Error: multiplying matrix A by state vector results in NaN matrix. \n");
		return;
	}


	mpfr_prec_t precRightSide;

	mpfr_t *Ax_sub_Bu;
	Ax_sub_Bu = allocateMPFRMatrix(1, n, 64);		//we set just 64 bits precision, but it will be changed in _exact functions

	mpfr_t *tmp2_n;
	tmp2_n = allocateMPFRMatrix(1, n, 64);

	MatrixSub_exact(Ax_sub_Bu, Ax, Bu, 1, n);
	MatrixSub_exact(tmp2_n, x_up, Ax_sub_Bu, 1, n);		//tmp2 = x_up - Ax - Bu

	for(i = 0; i < n; ++i)
	{		
		precRightSide = mpfr_get_prec(tmp2_n[i]);
		mpfr_set_prec(rightSide[i], precRightSide);
		mpfr_set_prec(rightSide[n + p + i], precRightSide);

		mpfr_set(rightSide[i], tmp2_n[i], MPFR_RNDN);
		mpfr_set(rightSide[n + p + i], tmp2_n[i], MPFR_RNDN);
		mpfr_neg(rightSide[n + p + i], rightSide[n + p + i], MPFR_RNDN);
	}

	mpfr_t *Cx_sub_Du;
	Cx_sub_Du = allocateMPFRMatrix(1, p, 64);

	mpfr_t *tmp2_p;
	tmp2_p = allocateMPFRMatrix(1, p, 64);

	MatrixSub_exact(Cx_sub_Du, Cx, Du, 1, p);
	MatrixSub_exact(tmp2_p, y_up, Cx_sub_Du, 1, p);		//tmp2 = x_up - Ax - Bu

	j = 0;
	for(i = n; i < n + p; ++i)
	{		
		precRightSide = mpfr_get_prec(tmp2_p[j]);
		mpfr_set_prec(rightSide[i], precRightSide);
		mpfr_set_prec(rightSide[n + p + i], precRightSide);
		
		mpfr_set(rightSide[i], tmp2_n[j], MPFR_RNDN);
		mpfr_set(rightSide[n + p + i], tmp2_n[j], MPFR_RNDN);
		mpfr_neg(rightSide[n + p + i], rightSide[n + p + i], MPFR_RNDN);
		j++;
	}

	///

	MatrixSub_exact(tmp2_n, x_low, Ax_sub_Bu, 1, n);

	j = 0;
	for(i = 2*(n+p); i < n + 2*(n+p); ++i)
	{		
		precRightSide = mpfr_get_prec(tmp2_n[j]);
		mpfr_set_prec(rightSide[i], precRightSide);
		mpfr_set_prec(rightSide[n + p + i], precRightSide);

		mpfr_set(rightSide[i], tmp2_n[j], MPFR_RNDN);
		mpfr_set(rightSide[n + p + i], tmp2_n[j], MPFR_RNDN);
		mpfr_neg(rightSide[n + p + i], rightSide[n + p + i], MPFR_RNDN);
		j++;
	}


	MatrixSub_exact(tmp2_p, y_low, Cx_sub_Du, 1, p);		//tmp2 = x_up - Ax - Bu

	j = 0;
	for(i = 3*n + 2*p; i < 3*(n+p); ++i)
	{		
		precRightSide = mpfr_get_prec(tmp2_p[j]);
		mpfr_set_prec(rightSide[i], precRightSide);
		mpfr_set_prec(rightSide[n + p + i], precRightSide);
		
		mpfr_set(rightSide[i], tmp2_n[j], MPFR_RNDN);
		mpfr_set(rightSide[n + p + i], tmp2_n[j], MPFR_RNDN);
		mpfr_neg(rightSide[n + p + i], rightSide[n + p + i], MPFR_RNDN);
		j++;
	}

	MatrixSub_exact(tmp2_n, x_up, x_low, 1, n);

	j = 0;
	for(i = 4*(n+p); i < 4*(n+p) + n; ++i)
	{
		precRightSide = mpfr_get_prec(tmp2_n[j]);
		mpfr_set_prec(rightSide[i], precRightSide);
		mpfr_set(rightSide[i], tmp2_n[j], MPFR_RNDN);

	}


	mpfr_t *tmp_q;
	tmp_q = allocateMPFRMatrix(1, q, 64);
	MatrixSub_exact(tmp_q, u_up, u_low, 1, q);

	j = 0;
	for(i = 4*(n+p) + n; i < 4*(n+p) + n + q; ++i)
	{
		precRightSide = mpfr_get_prec(tmp_q[j]);
		mpfr_set_prec(rightSide[i], precRightSide);
		mpfr_set(rightSide[i], tmp_q[j], MPFR_RNDN);
		j++;
	}

	freeMPFRMatrix(tmp_q, 1, q);
	freeMPFRMatrix(tmp2_n, 1, n);
	freeMPFRMatrix(tmp2_p, 1, p);
	freeMPFRMatrix(Ax, 1, n);
	freeMPFRMatrix(Bu, 1, n);
	freeMPFRMatrix(Cx, 1, p);
	freeMPFRMatrix(Du, 1, p);
	freeMPFRMatrix(Ax_sub_Bu, 1, n);
	freeMPFRMatrix(Cx_sub_Du, 1, p);
}



void BigMatrixConstruct(mpfr_t *Malade,\
	mpfr_t *A,mpfr_t *B, mpfr_t *C, mpfr_t *D, mpfr_t *rightSide,\
	uint64_t size, uint64_t cols, uint64_t n, uint64_t p, uint64_t q)
{

	// void minOnLine(mpfr_exp_t *v, mpfr_t *A, uint64_t n, uint64_t m)
	mpfr_exp_t *minA = (mpfr_exp_t*)malloc(n * sizeof(mpfr_exp_t));
	mpfr_exp_t *minC = (mpfr_exp_t*)malloc(p * sizeof(mpfr_exp_t));
	mpfr_exp_t *minB = (mpfr_exp_t*)malloc(n * sizeof(mpfr_exp_t));
	mpfr_exp_t *minD = (mpfr_exp_t*)malloc(p * sizeof(mpfr_exp_t));
	mpfr_exp_t *minExp = (mpfr_exp_t*)malloc(size * sizeof(mpfr_exp_t));

	minOnLine(minA, A, n, n);
	minOnLine(minB, B, n, q);
	minOnLine(minC, C, p, n);
	minOnLine(minD, D, p, q);

	int i,j;
	

	mpfr_t tmp;
	mpfr_init(tmp);
	mpfr_t one;
	mpfr_init(one);
	mpfr_set_ui(one, 1, MPFR_RNDN);

	int i_m, j_m;
	j_m = 0;
	i_m = 0;

	mpfr_exp_t scaling = 0;
	mpfr_prec_t prec;

	//Filling out rows 0 through n;
	for(i = 0; i < n; ++i)
	{
		scaling = min3(minA[i], minB[i], exponentMpfr(rightSide[i_m + i]));
		scaling = -scaling;
		minExp[i_m + i] = scaling;
		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);
		for(j = 0; j < n; ++j)
		{
			
			prec = 2*mpfr_get_prec(A[i * n + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], A[i * n + j], scaling, MPFR_RNDN);
		}
		j_m = n;
		for(j = 0; j < q; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = B[i * q + j] * 2^scaling
			prec = 2*mpfr_get_prec(B[i * q + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], B[i * q + j], scaling, MPFR_RNDN);
		}

		j_m = n + q;
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);
		mpfr_neg(Malade[(i_m + i)*cols + j_m + i], Malade[(i_m + i)*cols + j_m + i],MPFR_RNDN);
	
		j_m = 0;
	}

	//Filling out rows n through p ;

	i_m += n;
	j_m = 0;
	for(i = 0; i < p; ++i)
	{
		scaling = min3(minC[i], minD[i], exponentMpfr(rightSide[i_m + i]));
		scaling = -scaling;
		minExp[i_m + i] = scaling;
		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);
		for(j = 0; j < n; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = A[i * n + j] * 2^scaling
			prec = 2*mpfr_get_prec(C[i * n + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], C[i * n + j], scaling, MPFR_RNDN);
			
		}
		j_m = n;
		for(j = 0; j < q; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = B[i * q + j] * 2^scaling
			prec = 2*mpfr_get_prec(D[i * q + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], D[i * q + j], scaling, MPFR_RNDN);
			
		}

		j_m = n + q + n;
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);
		mpfr_neg(Malade[(i_m + i)*cols + j_m + i],Malade[(i_m + i)*cols + j_m + i],MPFR_RNDN);		//-I
	
		j_m = 0;
	}

	//-------
	//Filling out rows n + p through n + p + n;
	i_m = n + p;
	for(i = 0; i < n; ++i)
	{
		scaling = min3(minA[i], minB[i], exponentMpfr(rightSide[i_m + i]));
		scaling = -scaling;
		minExp[i_m + i] = scaling;
		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);
		for(j = 0; j < n; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = A[i * n + j] * 2^scaling
			prec = 2*mpfr_get_prec(A[i * n + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], A[i * n + j], scaling, MPFR_RNDN);
			mpfr_neg(Malade[(i_m + i)*cols + j_m + j],Malade[(i_m + i)*cols + j_m + j],MPFR_RNDN);		//-A
		}
		j_m = n;
		for(j = 0; j < q; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = B[i * q + j] * 2^scaling
			prec = 2*mpfr_get_prec(B[i * q + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], B[i * q + j], scaling, MPFR_RNDN);
			mpfr_neg(Malade[(i_m + i)*cols + j_m + j],Malade[(i_m + i)*cols + j_m + j],MPFR_RNDN);		//-B
		}

		j_m = n + q;
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);		//I
		
		j_m = 0;
	}

	//Filling out rows n + p + n through  p + n + p + n ;

	i_m += n;
	j_m = 0;
	for(i = 0; i < p; ++i)
	{
		scaling = min3(minC[i], minD[i], exponentMpfr(rightSide[i_m + i]));
		scaling = -scaling;
		minExp[i_m + i] = scaling;
		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);
		for(j = 0; j < n; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = A[i * n + j] * 2^scaling
			prec = 2*mpfr_get_prec(C[i * n + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], C[i * n + j], scaling, MPFR_RNDN);
			mpfr_neg(Malade[(i_m + i)*cols + j_m + j], Malade[(i_m + i)*cols + j_m + j],MPFR_RNDN);		//-C
		}
		j_m = n;
		for(j = 0; j < q; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = B[i * q + j] * 2^scaling
			prec = 2*mpfr_get_prec(D[i * q + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], D[i * q + j], scaling, MPFR_RNDN);
			mpfr_neg(Malade[(i_m + i)*cols + j_m + j], Malade[(i_m + i)*cols + j_m + j], MPFR_RNDN); //-D
		}

		j_m = n + q + n;
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);					//I

		j_m = 0;
	}

	//-----------
		//Filling out rows 2*n + 2*p  through 2*n + 2*p + n;
	i_m = 2*n + 2*p;
	for(i = 0; i < n; ++i)
	{
		scaling = min3(minA[i], minB[i], exponentMpfr(rightSide[i_m + i]));
		scaling = -scaling;
		minExp[i_m + i] = scaling;
		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);
		for(j = 0; j < n; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = A[i * n + j] * 2^scaling
			prec = 2*mpfr_get_prec(A[i * n + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], A[i * n + j], scaling, MPFR_RNDN);
		}
		j_m = n;
		for(j = 0; j < q; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = B[i * q + j] * 2^scaling
			prec = 2*mpfr_get_prec(B[i * q + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], B[i * q + j], scaling, MPFR_RNDN);
		}

		j_m = n + q + n + p;
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);
	
		j_m = 0;
	}

	//Filling out rows 2*n + 2*p + n through p + 2*n + 2*p + n ;

	i_m += n;
	j_m = 0;
	for(i = 0; i < p; ++i)
	{
		scaling = min3(minC[i], minD[i], exponentMpfr(rightSide[i_m + i]));
		scaling = -scaling;
		minExp[i_m + i] = scaling;
		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);

		for(j = 0; j < n; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = A[i * n + j] * 2^scaling
			prec = 2*mpfr_get_prec(C[i * n + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], C[i * n + j], scaling, MPFR_RNDN);
			
		}
		j_m = n;
		for(j = 0; j < q; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = B[i * q + j] * 2^scaling
			prec = 2*mpfr_get_prec(D[i * q + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], D[i * q + j], scaling, MPFR_RNDN);
		}

		j_m = n + q + n + p + n;
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);
	
		j_m = 0;
	}


	//---------

	//Filling out rows 3*n + 3*p  through 3*n + 3*p + n;
	i_m = 3*n + 3*p;
	for(i = 0; i < n; ++i)
	{
		scaling = min3(minA[i], minB[i], exponentMpfr(rightSide[i_m + i]));
		scaling = -scaling;
		minExp[i_m + i] = scaling;
		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);
		for(j = 0; j < n; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = A[i * n + j] * 2^scaling
			prec = 2*mpfr_get_prec(A[i * n + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], A[i * n + j], scaling, MPFR_RNDN);
			mpfr_neg(Malade[(i_m + i)*cols + j_m + j],Malade[(i_m + i)*cols + j_m + j],MPFR_RNDN);		//-A
		}
		j_m = n;
		for(j = 0; j < q; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = B[i * q + j] * 2^scaling
			prec = 2*mpfr_get_prec(B[i * q + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], B[i * q + j], scaling, MPFR_RNDN);
			mpfr_neg(Malade[(i_m + i)*cols + j_m + j],Malade[(i_m + i)*cols + j_m + j],MPFR_RNDN);		//-B
		}

		j_m = n + q + n + p;
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);		//-I
		mpfr_neg(Malade[(i_m + i)*cols + j_m + i],Malade[(i_m + i)*cols + j_m + i],MPFR_RNDN);
		j_m = 0;
	}

	//Filling out rows 3*n + 3*p + n through  3*n + 3*p + n + p ;

	i_m += n;
	j_m = 0;
	for(i = 0; i < p; ++i)
	{
		scaling = min3(minC[i], minD[i], exponentMpfr(rightSide[i_m + i]));
		scaling = -scaling;
		minExp[i_m + i] = scaling;

		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);

		for(j = 0; j < n; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = A[i * n + j] * 2^scaling
			prec = 2*mpfr_get_prec(C[i * n + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], C[i * n + j], scaling, MPFR_RNDN);
			mpfr_neg(Malade[(i_m + i)*cols + j_m + j], Malade[(i_m + i)*cols + j_m + j],MPFR_RNDN);		//-C
		}
		j_m = n;
		for(j = 0; j < q; ++j)
		{
			//Malade[(i_m + i)*cols + j_m + j = B[i * q + j] * 2^scaling
			prec = 2*mpfr_get_prec(D[i * q + j]);
			mpfr_set_prec(Malade[(i_m + i)*cols + j_m + j], prec);
			mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + j], D[i * q + j], scaling, MPFR_RNDN);
			mpfr_neg(Malade[(i_m + i)*cols + j_m + j], Malade[(i_m + i)*cols + j_m + j], MPFR_RNDN); //-D
		}

		j_m = n + q + n + p + n;
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);					//-I
		mpfr_neg(Malade[(i_m + i)*cols + j_m + i],Malade[(i_m + i)*cols + j_m + i],MPFR_RNDN);
		j_m = 0;
	}

	//--------
	//Filling out rows 4*n + 4*p through  4*n + 4*p + n ;
	i_m = 4*n + 4*p;
	j_m = 0;
	for(i = 0; i < n; ++i)
	{
		scaling = exponentMpfr(rightSide[i_m + i]);
		scaling = -scaling;
		minExp[i_m + i] = scaling;

		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);
	}
	i_m = 4*n + 4*p + n;
	j_m = n;
	for(i = 0; i < q; ++i)
	{
		scaling = exponentMpfr(rightSide[i_m + i]);
		scaling = -scaling;
		minExp[i_m + i] = scaling;

		mpfr_mul_2si(rightSide[i_m + i], rightSide[i_m + i],scaling, MPFR_RNDN);
		mpfr_mul_2si(Malade[(i_m + i)*cols + j_m + i], one, scaling, MPFR_RNDN);
	}


	// free(minA);
	// free(minB);
	// free(minC);
	// free(minD);
	// mpfr_clear(one);
	// mpfr_clear(tmp);
	// free(minExp);

}




void PrintLP(FILE *f, mpfr_t *A, \
			 mpfr_t *B,\
			 mpfr_t *C, \
			 mpfr_t *D, \
			 uint64_t n, uint64_t p, uint64_t q,\
			 mpfr_t *x_low, mpfr_t *x_up,\
			 mpfr_t *y_low, mpfr_t *y_up,\
			 mpfr_t *u_low, mpfr_t *u_up)
{


	mpfr_prec_t prec= 64;

 	uint64_t size = 5*n+4*p+q;
	uint64_t cols = 3 * n + 2 * p + q;
	mpfr_t *Malade = allocateMPFRMatrix(size, cols, 200);

	int i,j;
	for(i = 0; i < size; ++i)
	{
		for(j = 0; j < cols; ++j)
		{
			mpfr_set_ui(Malade[i * cols + j], (unsigned long int)0, MPFR_RNDN);
		}
	}


	mpfr_t *rightSide;
	rightSide = allocateMPFRMatrix(size, 1, 64);
	// printf("In function PrintLP, right side: \n");
	LPcomputeRightSide(rightSide, A, B, C, D, n, p, q,x_low, x_up,y_low, y_up, u_low, u_up);
	// writeMPFRMatrix(stderr, rightSide, size, 1, 100, MPFR_RNDN);

	BigMatrixConstruct(Malade, A,B, C, D, rightSide,size,cols, n, p, q);

	 // printf("\n \n \n Enormous matrix\n");
	 // writeMPFRMatrix(stderr, Malade, size, cols, 15, MPFR_RNDN);
	// printf("\n \n \n Right matrix\n");
	// writeMPFRMatrix(stderr, rightSide, size, 1, 100, MPFR_RNDN);
	


	fprintf(f, "Maximize\n");
	fprintf(f, " obj: ");

	char *strtmp = (char*)malloc(100000*sizeof(char));
	char *objstr= (char*)malloc(1000*sizeof(char));
	char *cstr= (char*)malloc(1000000*sizeof(char));


	for(i = 1; i < 2*n + 2*p; ++i)		//printing object function
	{
		sprintf(strtmp, "x%llu", i + n + q);
		strcat(objstr, strtmp);
		strcat(objstr, "+");
	}
	sprintf(strtmp, "x%llu", 2*n + 2*p);
	strcat(objstr, strtmp);


	fprintf(f, "%s\n", objstr);
	fprintf(f, "Subject To\n");
	mpfr_exp_t expptr;
	int nbits;

	mpz_t z;
	mpz_init2(z, 1000);

	int plus_flag = 0;

	for(i = 0; i < size; ++i)
	{
		fprintf(f, " c%d: ", i + 1);
		for(j = 0; j < cols; ++j)
		{
			if(mpfr_cmp_ui(Malade[i * cols + j], (unsigned long int)0) != 0)
			{
				mpfr_get_z(z, Malade[i * cols + j], MPFR_RNDN);
				// strtmp = mpz_get_str(strtmp, 10, z);
				// strcat(cstr, strtmp);
				// sprintf(strtmp, "x%d", j + 1);
				// strcat(cstr, strtmp);
				// strcat(cstr, "+");	
				mpz_out_str(f, 10, z);
				fprintf(f, " x%d+", j + 1);		
			}
		}

		fprintf(f, "0 < ");
		mpfr_get_z(z, rightSide[i], MPFR_RNDN);
		mpz_out_str(f, 10, z);
	
		fprintf(f, "\n");

	}

	fprintf(f, "General \n");
	for(i = 1; i <= cols; ++i)
		fprintf(f, "x%d\n",i);

	fprintf(f, "End\n");

}



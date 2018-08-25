/*
 ============================================================================
 Name        : matrixalg.c
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : This file contains glue and support functions for FXPF project.	  

 ============================================================================
*/

#include "matrixalg.h"

 

 void DoubleMatrixCopy(double *out, double *in, uint64_t n, uint64_t m)
 {
 	if(n == 0 || m == 0)
 	{
 		fprintf(stderr, "Error while copying matrices, incorrect size! \n");
 		return;
 	}
 	//int bytes = n * m * sizeof(double);
 	//memcpy(in, out, bytes); 		
 	int i,j;
 	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			out[i * m + j] = in[i * m + j];
		}
	}

 }

  void MPFRMatrixCopy(mpfr_t *out, mpfr_t *in, uint64_t n, uint64_t m)
  {
	if(n == 0 || m == 0)
 	{
 		fprintf(stderr, "Error while copying matrices, incorrect size! \n");
 		return;
 	}
 	mpfr_prec_t prec;
 	int i,j;
 	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			prec = mpfr_get_prec(in[i * m + j]) >= mpfr_get_prec(out[i * m + j]) ? mpfr_get_prec(in[i * m + j]) : mpfr_get_prec(out[i * m + j]);
			mpfr_set_prec(out[i * m + j], prec);
			mpfr_set(out[i * m + j], in[i * m + j], MPFR_RNDN);
		}
	}
  }

  void UINT64MatrixCopy(uint64_t *out, uint64_t *in, uint64_t n, uint64_t m)
 {
 	if(n == 0 || m == 0)
 	{
 		fprintf(stderr, "Error while copying matrices, incorrect size! \n");
 		return;
 	}	
 	int i,j;
 	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			out[i * m + j] = in[i * m + j];
		}
	}

 }

 int UINT64_to_INT(int *out, uint64_t *in, uint64_t n, uint64_t m)
 {
 	int i,j;
 	for(i = 0; i < n; ++i)
 	{
 		for(j = 0; j > m; ++j)
 		{
 			if(in[i * m + j] - INT_MIN <= (uint64_t)INT_MAX - INT_MIN)
 				out[i * m + j] = in[i * m + j];
 			else
 				fprintf(stderr, "Could not convert uint64 matrix to int \n");
 			return 0;
 		}
 	}
 	return 1;
 }

   void INT64MatrixCopy(int64_t *out, int64_t *in, uint64_t n, uint64_t m)
 {
 	if(n == 0 || m == 0)
 	{
 		fprintf(stderr, "Error while copying matrices, incorrect size! \n");
 		return;
 	}	
 	int i,j;
 	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			out[i * m + j] = in[i * m + j];
		}
	}

 }

 void INTMatrixCopy(int *out, int *in, uint64_t n, uint64_t m)
 {
 	 if(n == 0 || m == 0)
 	{
 		fprintf(stderr, "Error while copying matrices, incorrect size! \n");
 		return;
 	}	
 	int i,j;
 	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			out[i * m + j] = in[i * m + j];
		}
	}
 }

void DoubleMatrixPrint(FILE *file, double *A, uint64_t n, uint64_t m)
{
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			if(fprintf(file,"%f\t", A[i * m + j]) < 0)
			{
				fprintf(stderr, "Error writing matrix to a file\n");
				return;
			}
		}
		fprintf(file,"\n");
	}
}

void UINT64MatrixPrint(FILE *file, uint64_t *A, uint64_t n, uint64_t m)
{
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			if(fprintf(file,"%llu \t", A[i * m + j]) < 0)
			{
				fprintf(stderr, "Error writing matrix to a file\n");
				return;
			}
		}
		fprintf(file,"\n");
	}
}

  void INTMatrixPrint(FILE *file, int *A, uint64_t n, uint64_t m)
  {
  	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			if(fprintf(file,"%d \t", A[i * m + j]) < 0)
			{
				fprintf(stderr, "Error writing matrix to a file\n");
				return;
			}
		}
		fprintf(file,"\n");
	}
  }


void INT64MatrixPrint(FILE *file, int64_t *A, uint64_t n, uint64_t m)
{
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			if(fprintf(file,"%lld \t", A[i * m + j]) < 0)
			{
				fprintf(stderr, "Error writing matrix to a file\n");
				return;
			}
		}
		fprintf(file,"\n");
	}
}

void DoubleMatrixIdent(double *A, uint64_t n)
{
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			A[i * n + j] = 0.0;
		}
		A[i * n + i] = 1.0;
	}
}

void DoubleMatrixZero(double *A, uint64_t n, uint64_t m)
{
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			A[i * m + j] = 0.0;
		}
	}
}

void writeMPFRMatrix(FILE *stream, mpfr_t *A, uint64_t m, uint64_t n,size_t nmbr, mpfr_rnd_t rnd)
{
	fprintf(stream, "%llu %llu \n", m, n);
	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      mpfr_out_str(stream, (int)10, nmbr, A[i * n + j], rnd);
		      fprintf(stream, "\t");
	      }
	      fprintf(stream, "\n");
	}

}

void writeMPFIMatrix(FILE *stream, mpfi_t *A, uint64_t m, uint64_t n,size_t nmbr)
{
	fprintf(stream, "%llu %llu \n", m, n);
	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      mpfi_out_str(stream, (int)10, nmbr, A[i * n + j]);
		      fprintf(stream, "\t");
	      }
	      fprintf(stream, "\n");
	}

}


void readMPFRMatrix(mpfr_t *A, FILE *stream, uint64_t m, uint64_t n, mpfr_rnd_t rnd)
{

	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      mpfr_inp_str(A[i * n + j], stream, (int)10, rnd);
	      }
	}
}

void readDoubleMatrix(FILE *stream, double *A,uint64_t m, uint64_t n)
{
	
	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      fscanf(stream, " %lf ", &(A[i*n + j]));			
	      }
	}
  
}

void readINTMatrix(FILE *stream, int *A,uint64_t m, uint64_t n)
{
	
	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      fscanf(stream, " %d ", &(A[i*n + j]));			
	      }
	}
  
}

void doubleToMPFRMatrix(mpfr_t *mpA, double *A, uint64_t m, uint64_t n)
{
	int i, j;
	mp_prec_t prec = 70;

	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			mpfr_set_prec(mpA[i * n + j], prec);
			mpfr_set_d(mpA[i * n + j], A[i * n + j], MPFR_RNDN);
		}
	}
}

/* Allocates a n * m matrix and initializes all entries to prec
   bits 
*/
mpfr_t *allocateMPFRMatrix(uint64_t n,uint64_t m, mp_prec_t prec) {
  mpfr_t *A;
  uint64_t i, j;

  A = (mpfr_t *) calloc(n * m, sizeof(mpfr_t));
  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfr_init2(A[i * m + j], prec);
    }
  }

  return A;
}



/* Clears all entries of a n * m matrix A and frees the memory
   allocated to that matrix 
*/
void freeMPFRMatrix(mpfr_t *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfr_clear(A[i * m + j]);
    }
  }
 
  free(A);
}


/* Allocates a n sized vector and initializes all entries to prec
   bits 
*/
mpfr_t *allocateMPFRVector(uint64_t n, mp_prec_t prec) {
  mpfr_t *v;
  uint64_t i;

  v = (mpfr_t *) calloc(n, sizeof(mpfr_t));
  for (i=0;i<n;i++) {
    mpfr_init2(v[i], prec);
  }

  return v;
}

/* Clears all entries of a n sized vector v and frees the memory
   allocated to that vector
*/
void freeMPFRVector(mpfr_t *v, uint64_t n) {
  uint64_t i;

  for (i=0;i<n;i++) {
      mpfr_clear(v[i]);
  }
 
  free(v);
}

/* Sets a n * m matrix A to all zeros */
void setMatrixZero(mpfr_t *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfr_set_si(A[i * m + j], 0, GMP_RNDN); /* exact */
    }
  }
}

/* Sets a n * n matrix A to identity */
void setMatrixIdentity(mpfr_t *A, uint64_t n) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      mpfr_set_si(A[i * n + j], 0, GMP_RNDN); /* exact */
    }
    mpfr_set_si(A[i * n + i], 1, GMP_RNDN); /* exact */
  }
}

/* Vertically concatanates two double matrices.
Input: A is m x n matrix, B is k x n matrix.
Output: C is (m + k) x n matrix.
Matrix C is assumed to be preallocated outside the function. */
void doubleMatrixConcatVertical(double *C, double *A, double *B, uint64_t m, uint64_t k, uint64_t n)
{
	int i,j;
	for(j = 0; j < n; ++j)
	{
		for(i = 0; i < m; ++i)
		{
			C[i * n + j] = A[i * n + j];
		}
		for(i = 0; i < k; ++i)
		{
			C[(i + m) * n + j] = B[i * n + j];
		}
	}
}

/* Horizontally concatanates two double matrices.
Input: A is m x n matrix, B is m x k matrix.
Output: C is m x (n + k) matrix.
Matrix C is assumed to be preallocated outside the function. */
void doubleMatrixConcatHorizontal(double *C, double *A, double *B, uint64_t m, uint64_t n, uint64_t k)
{
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			C[i * (n + k) + j] = A[i * n + j];
		}
		for(j = 0; j < k; ++j)
		{
			C[i * (n + k) + (j + n)] = B[i * k + j];
		}
	}
}

mp_prec_t getMaxPrec(mpfr_t op1, mpfr_t op2)
{
	mp_prec_t prec1 = mpfr_get_prec(op1);
	mp_prec_t prec2 = mpfr_get_prec(op2);
	return prec1 >= prec2 ? prec1 : prec2;
}

/* Horizontally concatanates two MPFR matrices.
Input: A is m x n matrix, B is m x k matrix.
Output: C is m x (n + k) matrix.
Matrix C is assumed to be preallocated outside the function. */
void MPFRMatrixConcatHorizontal(mpfr_t *C, mpfr_t *A, mpfr_t *B, uint64_t m, uint64_t n, uint64_t k)
{
	int i,j;
	mp_prec_t maxPrec;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			maxPrec = getMaxPrec(A[i * n + j], C[i * (n + k) + j]);
			mpfr_set_prec(C[i * (n + k) + j], maxPrec);
			mpfr_set(C[i * (n + k) + j],A[i * n + j], MPFR_RNDN);				//C[i * (n + k) + j] = A[i * n + j];
		}
		for(j = 0; j < k; ++j)
		{
			maxPrec = getMaxPrec(B[i * n + j], C[i * (n + k) + j]);
			mpfr_set_prec(C[i * (n + k) + (j + n)], maxPrec);
			mpfr_set(C[i * (n + k) + (j + n)],B[i * k + j], MPFR_RNDN);			//C[i * (n + k) + (j + n)] = B[i * k + j];
		}
	}
}

/* Vertically concatanates two MPFR matrices.
Input: A is m x n matrix, B is k x n matrix.
Output: C is (m + k) x n matrix.
Matrix C is assumed to be preallocated outside the function. */
void MPFRMatrixConcatVertical(mpfr_t *C, mpfr_t *A, mpfr_t *B, uint64_t m, uint64_t k, uint64_t n)
{
	int i,j;
	mp_prec_t maxPrec;
	for(j = 0; j < n; ++j)
	{
		for(i = 0; i < m; ++i)
		{
			maxPrec = getMaxPrec(A[i * n + j], C[i * n + j]);
			mpfr_set_prec(C[i * n + j], maxPrec);
			mpfr_set(C[i * n + j],A[i * n + j], MPFR_RNDN); 					//C[i * n + j] = A[i * n + j];
		}
		for(i = 0; i < k; ++i)
		{
			maxPrec = getMaxPrec(C[(i + m) * n + j], B[i * n + j]);
			mpfr_set_prec(C[(i + m) * n + j], maxPrec);
			mpfr_set(C[(i + m) * n + j], B[i * n + j], MPFR_RNDN); 				//C[(i + m) * n + j] = B[i * n + j];
		}
	}
}

void ConcatTest(double *A, uint64_t n)
{
	printf("testing Matrix vertical concat\n");
	double  *Z;
	Z = (double*)malloc(4 * n * sizeof(double));
	DoubleMatrixZero(Z, 4, n);
	double  *Concat;
	Concat = (double*)malloc((4 + n) * n * sizeof(double));

	
	doubleMatrixConcatVertical(Concat, A, Z, n, 4, n);
	DoubleMatrixPrint(stderr, Concat, n+4, n);

	printf("testing Matrix horisontal concat\n");
	double  *I;
	I = (double*)malloc(n * n * sizeof(double));
	DoubleMatrixIdent(I, n);
	double  *ConcatH;
	ConcatH = (double*)malloc((n + n) * n * sizeof(double));

	
	doubleMatrixConcatHorizontal(ConcatH, A, I, n, n, n);
	DoubleMatrixPrint(stderr, ConcatH, n, n + n);
}

mp_prec_t getMaxPrecMPFRMatrix(mpfr_t *A, uint64_t m, uint64_t n)
{
	int i,j;
	mp_prec_t current = 0;
	mp_prec_t max = 0;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			current = mpfr_get_prec(A[i * n + j]);
			if(current >= max) max = current;
		}
	}
	return max;
}

/* CHeck if any of the elements of a double m x n matrix is NaN.
Returns a non-zero value (true) if A has NaN value, and zero (false) otherwise. */
int matrixIsNan_double(double *A, uint64_t m, uint64_t n)
{
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
		    if(isnan(A[i*n + j])) return 1;	
		}
	}
	return 0;
}

void setAllto(int* A, int val, uint64_t m, uint64_t n)
{
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			A[i * n + j] = val;
		}
	}

}


/* CHeck if any of the elements of a double m x n matrix is NaN.
Returns a non-zero value (true) if A has NaN value, and zero (false) otherwise. */
int matrixIsNan_mpfr(mpfr_t *A, uint64_t m, uint64_t n)
{
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
		    if(!mpfr_number_p(A[i*n + j])) return 1;	
		}
	}
	return 0;

}

void mpfrToUINT64(uint64_t *B, mpfr_t *A, uint64_t m, uint64_t n, mpfr_rnd_t rnd)
{
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
		    B[i * n + j] = mpfr_get_ui(A[i*n + j], rnd);	
		}
	}
}

mpfr_exp_t min3(mpfr_exp_t e1, mpfr_exp_t e2, mpfr_exp_t e3)
{
	mpfr_exp_t res = e1 < e2 ? e1 : e2;
	return res < e3 ? res : e3;
}



/* Returns the exact product of op1 and op2 in the variable res.
Changes the precision of res. */
void mul_exact(mpfr_t res, mpfr_t op1, mpfr_t op2)
{
	if(!mpfr_number_p(op1) || !mpfr_number_p(op2))
	{
		fprintf(stderr, "Trying to multiply two numbers, one of which is NaN(Inf)\n");
		mpfr_set_zero(res, 1);
		return;
	}
	mpfr_prec_t res_prec = mpfr_get_prec(op1) + mpfr_get_prec(op2) + 5;
	mpfr_set_prec(res, res_prec);
	mpfr_mul(res, op1, op2, MPFR_RNDN);
}

/* Returns the exact product of op1 and op2 in the variable res.
Changes the precision of res  */
void sum_exact(mpfr_t res, mpfr_t op1, mpfr_t op2)
{
	if(!mpfr_number_p(op1) || !mpfr_number_p(op2))
	{
		fprintf(stderr, "Trying to add two numbers, one of which is NaN(Inf). \n");
		mpfr_set_zero(res, 1);
		return;
	}
	mpfr_prec_t op_prec_max = mpfr_get_prec(op1) >= mpfr_get_prec(op2) ? mpfr_get_prec(op1) : mpfr_get_prec(op2);
	mpfr_prec_t res_prec = op_prec_max + ceil(log2(op_prec_max)) + 5;
	mpfr_set_prec(res, res_prec);
	if(mpfr_add(res, op1, op2, MPFR_RNDN) != 0)
		fprintf(stderr, "Warning: exact summation is not exact..\n");
}


/*
	The function performs an exact multiplication of an MPFR matrix by a vector.
	Matrix is of size n x m, vector is of size m x 1.
	The space for the result is assumed preallocated.
	The precision of the result is changed in order to hold the exact sum.
 */
void MatrixMulVector_exact(mpfr_t *res, mpfr_t *A, mpfr_t *x, uint64_t n, uint64_t m, mpfr_t *scratch2)
{

	/*
	The srategy is simple. Denote B: = A * x
		B_i = sum_{k=0}^m A_ik * x_k

	We compute each product A_ik * x_k exactly and store it in variable t1.
	The precision of t1 is changed each time before multiplication inside the
	mul_exact function.

	Then, we sum the product t1 exactly with B_i and store the result in t2.
	The precision of t2 is changed each time before summaion inside the
	sum_exact function.

	Finally, we adjust the precision of B_i such that it can hold exactly the
	value from t2 and copy the value to the B_i variable.
	 */

/*	mpfr_t t1, t2;
	mpfr_init(t1);
	mpfr_init(t2);*/


	mpfr_prec_t t2_prec;
	int i,k;
	for(i = 0; i < n; ++i)
	{
		mpfr_set_ui(res[i], 0, MPFR_RNDN);
		for(k=0; k < m; ++k)
		{
			mul_exact(scratch2[0], A[i * m + k], x[k]);		// t1 = A_ik*x_k exact
			sum_exact(scratch2[1], scratch2[0], res[i]);				// t2 = t1 + res[i] exact
			if(!mpfr_number_p(scratch2[1]))
				fprintf(stderr, "Trying to multiply two matrices, which have a NaN/Inf element\n");
			t2_prec = mpfr_get_prec(scratch2[1])+1;			// set precision of res[i] to precision of t2
			mpfr_set_prec(res[i], t2_prec);
			mpfr_set(res[i], scratch2[1], MPFR_RNDN);		//res = res + A_ik*x_k exact
			
		}
	}
	/*mpfr_clear(t1);
	mpfr_clear(t2);*/
}

/* res = A[line, :] */
void MatrixExtractRow(mpfr_t *res, mpfr_t *A, uint64_t n, uint64_t m, uint64_t row)
{
	
	int j;
	mpfr_prec_t prec;
	for(j = 0; j < m; ++j)
	{
		prec = mpfr_get_prec(A[row * m + j]);
		mpfr_set_prec(res[j], prec);
		mpfr_set(res[j], A[row * m + j], MPFR_RNDN);
	}

}

/* res = A[:, col] */
void MatrixExtractColumn(mpfr_t *res, mpfr_t *A, uint64_t n, uint64_t m, uint64_t col)
{
	int i;
	mpfr_prec_t prec;
	for(i = 0; i < n; ++i)
	{
		prec = mpfr_get_prec(A[i * m + col]);
		mpfr_set_prec(res[i], prec);
		mpfr_set(res[i], A[i * m + col], MPFR_RNDN);
	
	}

}

/* For a T x n matrix A sets its k-th row values to the values of
1 x n vector v.
Matrix A is assumed preallocated. The precision of its elements may be changed 
in order to containt the values of vector v.
  */
void MatrixSetRow(mpfr_t *A, mpfr_t *v, uint64_t T, uint64_t n, uint64_t k)
{
	int i,j;
	mpfr_prec_t prec;
	for(i = 0; i < n; ++i)
	{
		/* A[k * n + i]  = v[i]*/
		prec = mpfr_get_prec(A[k * n + i]) >= mpfr_get_prec(v[i]) ? mpfr_get_prec(A[k * n + i]) : mpfr_get_prec(v[i]);
		mpfr_set_prec(A[k * n + i], prec);
		mpfr_set(A[k * n + i], v[i], MPFR_RNDN);
	}
}

/*
	The function performs an exact summation of two MPFR matrices.
	Each matrix is of size n x m.
	The space for the result is assumed preallocated.
	The precision of the result is changed in order to hold the exact sum.
	The pointer to res must be distinct of pointer of other arguments. 
 */
void MatrixAdd_exact(mpfr_t *res, mpfr_t *A, mpfr_t *B, uint64_t n, uint64_t m)
{
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			sum_exact(res[i * m + j], A[i * m + j], B[i * m + j]);
		}
	}
}

void MatrixSub_exact(mpfr_t *res, mpfr_t *A, mpfr_t *B, uint64_t n, uint64_t m)
{
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			mpfr_neg(B[i * m + j], B[i * m + j], MPFR_RNDN);
			sum_exact(res[i * m + j], A[i * m + j], B[i * m + j]);
		}
	}
}

void MatrixQuantize(mpfr_t *res, mpfr_t *v, int *lsb, uint64_t n, uint64_t m, mpfr_t q)
{
	int i, j;
	//mpfr_t q;
	//mpfr_init(q);
	mpfr_prec_t prec;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{

			// printf("i = %d j = %d\n",i, j );
			prec = mpfr_get_prec(v[i * m + j]);
			mpfr_set_prec(q, prec);
			mpfr_set_prec(res[i * m + j], prec);
			/* q = v[i] * 2^(-l) */
			/* usually l < 0, however it may be positive */

			/* if l <= 0, then q = v[i] / 2^l */
			if(lsb[i * m + j] <= 0)
			{
				mpfr_set(q, v[i * m + j], MPFR_RNDN);
				if(mpfr_mul_2ui(q, q, abs(lsb[i * m + j]), MPFR_RNDN) != 0)
				{
					fprintf(stderr, "Warning: in VectorQuantization one of the components is not correctly rounded \n");
				}

			}
			else
			{
				if(mpfr_div_2ui(q, v[i * m + j], abs(lsb[i * m + j]), MPFR_RNDN) != 0)
				{
					fprintf(stderr, "Warning: in VectorQuantization one of the components is not correctly rounded \n");
				}
			}

			/* q = round(q) */

			//mpfr_round(q, q);
			mpfr_rint(q, q, MPFR_RNDN);
			// mpfr_set(resINT[i * m + j], q, MPFR_RNDN);
		
			/* res[i] = q * 2^l*/
			if(lsb[i * m + j] <= 0)
			{

				if(mpfr_div_2ui(res[i * m + j], q, abs(lsb[i * m + j]), MPFR_RNDN) != 0)
				{
					fprintf(stderr, "Warning: in VectorQuantization one of the components is not correctly rounded \n");
				}

			}
			else
			{
				if(mpfr_mul_2ui(res[i * m + j], v[i * m + j], abs(lsb[i * m + j]), MPFR_RNDN) != 0)
				{
					fprintf(stderr, "Warning: in VectorQuantization one of the components is not correctly rounded \n");
				}
			}
		}
		

	}
	//mpfr_clear(q);
}

void VectorQuantize(mpfr_t *res, mpfr_t *v, int *lsb, uint64_t n, mpfr_t q)
{
	int i;
	//mpfr_t q;
	//mpfr_init(q);
	mpfr_prec_t prec;


	for(i = 0; i < n; ++i)
	{
		prec = mpfr_get_prec(v[i]);
		mpfr_set_prec(q, prec);
		mpfr_set_prec(res[i], prec);
		/* q = v[i] * 2^(-l) */
		/* usually l < 0, however it may be positive */

		/* if l <= 0, then q = v[i] * 2^-l =  v[i] * 2^abs(l)*/
		if(lsb[i] <= 0)
		{
			mpfr_set(q, v[i], MPFR_RNDD);
			if(mpfr_mul_2ui(q, q, abs(lsb[i]), MPFR_RNDN) != 0)
			{
				fprintf(stderr, "Warning: in VectorQuantization one of the components is not correctly rounded \n");
			}


		}
		else
		/* if l > 0, then q = v[i] * 2^-l =  v[i] / 2^abs(l)*/
		{
			if(mpfr_div_2ui(q, v[i], abs(lsb[i]), MPFR_RNDN) != 0)
			{
				fprintf(stderr, "Warning: in VectorQuantization one of the components is not correctly rounded \n");
			}
		}

		/* q = round(q) */

		//mpfr_round(q, q);
		mpfr_rint(q, q, MPFR_RNDD);
	
		/* res[i] = q * 2^l*/
		if(lsb[i] <= 0)
		{
			if(mpfr_div_2ui(res[i], q, abs(lsb[i]), MPFR_RNDN) != 0)
			{
				fprintf(stderr, "Warning: in VectorQuantization one of the components is not correctly rounded \n");
			}

		}
		else
		{
			if(mpfr_mul_2ui(res[i], v[i], abs(lsb[i]), MPFR_RNDN) != 0)
			{
				fprintf(stderr, "Warning: in VectorQuantization one of the components is not correctly rounded \n");
			}
		}

	}
	//mpfr_clear(q);
}

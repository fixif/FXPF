/*
 ============================================================================
 Name        : matrixalg.h
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : This file contains glue and support functions for FXPF project.	  

 ============================================================================
*/

#ifndef _MATRIXALG_H
#define _MATRIXALG_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <mpfi_io.h>
#include <limits.h>


/* Copies content of a nxm matrix B into matrix A of the same size. */
 void MPFRMatrixCopy(mpfr_t *out, mpfr_t *in, uint64_t n, uint64_t m);
void DoubleMatrixCopy(double *out, double *in, uint64_t n, uint64_t m);
  void INT64MatrixCopy(int64_t *out, int64_t *in, uint64_t n, uint64_t m);
    void INTMatrixCopy(int *out, int *in, uint64_t n, uint64_t m);

  void UINT64MatrixCopy(uint64_t *out, uint64_t *in, uint64_t n, uint64_t m);
  void INTMatrixPrint(FILE *file, int *A, uint64_t n, uint64_t m);
  void INT64MatrixPrint(FILE *file, int64_t *A, uint64_t n, uint64_t m);
/*Prints a n x m matrix to an output file. */
void DoubleMatrixPrint(FILE *file, double *A, uint64_t n, uint64_t m);
void UINT64MatrixPrint(FILE *file, uint64_t *A, uint64_t n, uint64_t m);

/* Sets a square matrix A to Identity matrix */
void DoubleMatrixIdent(double *A, uint64_t n);

int UINT64_to_INT(int *out, uint64_t *in, uint64_t n, uint64_t m);

/* Sets a n x m matrix A to zero matrix */
void DoubleMatrixZero(double *A, uint64_t n, uint64_t m);

/* Sets a n * m matrix A to all zeros */
void setMatrixZero(mpfr_t *A, uint64_t n, uint64_t m);

/* Sets a n * n matrix A to identity */
void setMatrixIdentity(mpfr_t *A, uint64_t n);

/*
Write to file stream a complex m * n matrix rounded in the direction rnd with its real and imaginary parts in ReA and ImA respectively.
The function prints nmbr significant digits exactly, or if nmbr is 0, enough digits
so that matrix could be read back exactly.
Format of output: first line is two digits, representing size of matrix.
Then values are printed separated with tabulation.
The function prints matrix in base 10.
*/
void writeMPFRMatrix(FILE *stream, mpfr_t *A, uint64_t m, uint64_t n,size_t nmbr, mpfr_rnd_t rnd);

/*
Write to file stream a complex m * n matrix rounded in the direction rnd with its real and imaginary parts in ReA and ImA respectively.
The function prints nmbr significant digits exactly, or if nmbr is 0, enough digits
so that matrix could be read back exactly.
Format of output: first line is two digits, representing size of matrix.
Then values are printed separated with tabulation.
The function prints matrix in base 10.
*/
void writeMPFIMatrix(FILE *stream, mpfi_t *A, uint64_t m, uint64_t n,size_t nmbr);
/*
Read a floating-point matrix A of size m * n from file stream, using rounding direction rnd.
Matrix A is assumed to be declared and initialized outside this function. Precision of matrix A is set outside this function.
Format of input: floating-point numbers must be in base 10 in form A@B or AeB, where A is mantissa and B is exponent.
*/
void readMPFRMatrix(mpfr_t *A, FILE *stream, uint64_t m, uint64_t n, mpfr_rnd_t rnd);

void readDoubleMatrix(FILE *stream, double *A,uint64_t m, uint64_t n);
void readINTMatrix(FILE *stream, int *A,uint64_t m, uint64_t n);

/* Convert a double n * m matrix A  to the MPFR matrix mpA.
Matrix ReA is assumed to be declared and pre-allocated outside the function.
THe function changes precision of ReA to 64. 
*/
void doubleToMPFRMatrix(mpfr_t *mpA, double *A, uint64_t m, uint64_t n);


/* Allocates a n * m matrix and initializes all entries to prec
   bits 
*/
mpfr_t *allocateMPFRMatrix(uint64_t n,uint64_t m, mp_prec_t prec);

/* Clears all entries of a n * m matrix A and frees the memory
   allocated to that matrix 
*/
void freeMPFRMatrix(mpfr_t *A, uint64_t n, uint64_t m);

/* Allocates a n sized vector and initializes all entries to prec
   bits 
*/
mpfr_t *allocateMPFRVector(uint64_t n, mp_prec_t prec);

/* Clears all entries of a n sized vector v and frees the memory
   allocated to that vector
*/
void freeMPFRVector(mpfr_t *v, uint64_t n);

/* Vertically concatanates two double matrices.
Input: A is m x n matrix, B is k x n matrix.
Output: C is (m + k) x n matrix.
Matrix C is assumed to be preallocated outside the function. */
void doubleMatrixConcatVertical(double *C, double *A, double *B, uint64_t m, uint64_t k, uint64_t n);


/* Vertically concatanates two MPFR matrices.
Input: A is m x n matrix, B is k x n matrix.
Output: C is (m + k) x n matrix.
Matrix C is assumed to be preallocated outside the function. */
void MPFRMatrixConcatVertical(mpfr_t *C, mpfr_t *A, mpfr_t *B, uint64_t m, uint64_t k, uint64_t n);


/* Horizontally concatanates two double matrices.
Input: A is m x n matrix, B is m x k matrix.
Output: C is m x (n + k) matrix.
Matrix C is assumed to be preallocated outside the function. */
void doubleMatrixConcatHorizontal(double *C, double *A, double *B, uint64_t m, uint64_t n, uint64_t k);

/* Horizontally concatanates two MPFR matrices.
Input: A is m x n matrix, B is m x k matrix.
Output: C is m x (n + k) matrix.
Matrix C is assumed to be preallocated outside the function. */
void MPFRMatrixConcatHorizontal(mpfr_t *C, mpfr_t *A, mpfr_t *B, uint64_t m, uint64_t n, uint64_t k);


/*Provides test for concat functions. To be removed. */
void ConcatTest(double *A, uint64_t n);

/* Returns max precision of two elements. */
mp_prec_t getMaxPrec(mpfr_t op1, mpfr_t op2);

mp_prec_t getMaxPrecMPFRMatrix(mpfr_t *A, uint64_t m, uint64_t n);

/* CHeck if any of the elements of a double m x n matrix is NaN.
Returns a non-zero value (true) if A has NaN value, and zero (false) otherwise. */
int matrixIsNan_double(double *A, uint64_t m, uint64_t n);

/* CHeck if any of the elements of a double m x n matrix is NaN.
Returns a non-zero value (true) if A has NaN value, and zero (false) otherwise. */
int matrixIsNan_mpfr(mpfr_t *A, uint64_t m, uint64_t n);

void mpfrToUINT64(uint64_t *B, mpfr_t *A, uint64_t m, uint64_t n, mpfr_rnd_t rnd);
void setAllto(int* A, int val, uint64_t m, uint64_t n);

mpfr_exp_t min3(mpfr_exp_t e1, mpfr_exp_t e2, mpfr_exp_t e3);




void mul_exact					(mpfr_t res, mpfr_t op1, mpfr_t op2);
void sum_exact					(mpfr_t res, mpfr_t op1, mpfr_t op2);
//void MatrixMulVector_exact		(mpfr_t *res, mpfr_t *A, mpfr_t *x, uint64_t n, uint64_t m);
void MatrixExtractRow			(mpfr_t *res, mpfr_t *A, uint64_t n, uint64_t m, uint64_t row);
void MatrixExtractColumn		(mpfr_t *res, mpfr_t *A, uint64_t n, uint64_t m, uint64_t col);
void MatrixSetRow				(mpfr_t *A, mpfr_t *v, uint64_t T, uint64_t n, uint64_t k);
void MatrixAdd_exact			(mpfr_t *res, mpfr_t *A, mpfr_t *B, uint64_t n, uint64_t m);
void MatrixSub_exact			(mpfr_t *res, mpfr_t *A, mpfr_t *B, uint64_t n, uint64_t m);
//void MatrixQuantize				(mpfr_t *res, mpfr_t *v, int *lsb, uint64_t n, uint64_t m);
//void VectorQuantize				(mpfr_t *res, mpfr_t *v, int *lsb, uint64_t n);

void VectorQuantize(mpfr_t *res, mpfr_t *v, int *lsb, uint64_t n, mpfr_t q);
void MatrixQuantize(mpfr_t *res, mpfr_t *v, int *lsb, uint64_t n, uint64_t m, mpfr_t q);
void MatrixMulVector_exact(mpfr_t *res, mpfr_t *A, mpfr_t *x, uint64_t n, uint64_t m, mpfr_t *scratch2);

#endif

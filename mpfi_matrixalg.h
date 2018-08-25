/*
 ============================================================================
 Name        : mpfi_matrixalg.h
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : This header file contains  prototypes for support functions 
 				for matrix representation in MPFI	  

 ============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpfr.h>
#include <mpfi.h>
#include <mpfi_io.h>

void *safeMalloc(size_t size);
void *safeCalloc(size_t nmemb, size_t size);
void *safeRealloc(void *ptr, size_t size);
void safeFree(void *ptr);

 /* Allocates a n * m matrix and initializes all entries to prec
   bits 
*/
mpfi_t *allocateMPFIMatrix(uint64_t n,uint64_t m, mp_prec_t prec);

/* Clears all entries of a n * m matrix A and frees the memory
   allocated to that matrix 
*/
void freeMPFIMatrix(mpfi_t *A, uint64_t n, uint64_t m);

/* Sets a n * m matrix A to all zeros */
void setMPFIMatrixZero(mpfi_t *A, uint64_t n, uint64_t m);

/* Sets a n * n matrix A to identity */
void setMPFIMatrixIdentity(mpfi_t *A, uint64_t n);

/* Converts element-by-element MPFR matrix A
into MPFI latrix B.
Matrix B is assumed preallocated outside the function. */
void MPFItoMPFRMatrix(mpfi_t *B, mpfr_t *A, uint64_t m, uint64_t n);

mp_prec_t getMaxPrecMPFIMatrix(mpfi_t *A, uint64_t m, uint64_t n);


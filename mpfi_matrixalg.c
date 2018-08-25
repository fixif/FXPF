/*
 ============================================================================
 Name        : mpfi_matrixalg.c
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : This file contains support functions for matrix representation
 				in MPFI	  

 ============================================================================
*/

 #include "mpfi_matrixalg.h"

/* Allocates size bytes of memory; on failure prints an error message
   and calls exit(1) 
*/
void *safeMalloc(size_t size) {
  void *ptr;
  
  ptr = malloc(size);
  if (ptr == NULL) {
    fprintf(stderr, "Could not allocate %zd bytes of memory\n", size);
    exit(1);
  }

  return ptr;
}

/* Allocates an array of nmemb elements, each of which occupies size
   bytes of memory; on failure prints an error message and calls
   exit(1)
*/
void *safeCalloc(size_t nmemb, size_t size) 
{
  void *ptr;
  
  ptr = calloc(nmemb, size);
  if (ptr == NULL) {
    fprintf(stderr, "Could not allocate an array of %zd elements with %zd bytes of memory each\n", nmemb, size);
    exit(1);
  }

  return ptr;
}

/* Reallocates size bytes of memory at the location pointed by ptr; 
   on failure prints an error message and calls exit(1)
*/
void *safeRealloc(void *ptr, size_t size) 
{
  void *newPtr;
  
  newPtr = realloc(ptr, size);
  if ((newPtr == NULL) && (!((size == 0) && (ptr != NULL)))) {
      fprintf(stderr, "Could not rellocate %zd bytes of memory at address %p\n", size, ptr);
    exit(1);
  }

  return newPtr;
}

/* Frees the memory at the location pointed by ptr */
void safeFree(void *ptr) {
  free(ptr);
}

 /* Allocates a n * m matrix and initializes all entries to prec
   bits 
*/
mpfi_t *allocateMPFIMatrix(uint64_t n,uint64_t m, mp_prec_t prec) {
  mpfi_t *A;
  uint64_t i, j;

  A = (mpfi_t *) safeCalloc(n * m, sizeof(mpfi_t));
  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfi_init2(A[i * m + j], prec);
    }
  }

  return A;
}


/* Clears all entries of a n * m matrix A and frees the memory
   allocated to that matrix 
*/
void freeMPFIMatrix(mpfi_t *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfi_clear(A[i * m + j]);
    }
  }
 
  safeFree(A);
}


/* Sets a n * m matrix A to all zeros */
void setMPFIMatrixZero(mpfi_t *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfi_set_si(A[i * m + j], 0); 
    }
  }
}

/* Sets a n * n matrix A to identity */
void setMPFIMatrixIdentity(mpfi_t *A, uint64_t n) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      mpfi_set_si(A[i * n + j], 0); 
    }
    mpfi_set_si(A[i * n + i], 1); 
  }
}

/* Converts element-by-element MPFR matrix A
into MPFI latrix B.
Matrix B is assumed preallocated outside the function. */
void MPFItoMPFRMatrix(mpfi_t *B, mpfr_t *A, uint64_t m, uint64_t n)
{
	int i,j;
	mp_prec_t prec;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			prec = mpfr_get_prec(A[i * n + j]);
			mpfi_set_prec(B[i * n + j], prec);
			mpfi_set_fr(B[i * n + j], A[i * n + j]);
		}
	}

}

mp_prec_t getMaxPrecMPFIMatrix(mpfi_t *A, uint64_t m, uint64_t n)
{
	int i,j;
	mp_prec_t current = 0;
	mp_prec_t max = 0;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			current = mpfi_get_prec(A[i * n + j]);
			if(current >= max) max = current;
		}
	}
	return max;
}







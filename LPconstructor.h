/*
 ============================================================================
 Name        : LPconstructor.h
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : This file contains header files for the module responsible 
 				for the consutruction of the Integer Linear Programming problem
 				and writing it into the *.lp format. (which is later to be used 
 				with SCIP)

 ============================================================================
*/

#ifndef _LPCONSTRUCTOR_H
#define _LPCONSTRUCTOR_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "matrixalg.h"

mp_exp_t exponentMpfr(mpfr_t x);

void PrintLP(FILE *f, mpfr_t *A, \
			 mpfr_t *B,\
			 mpfr_t *C, \
			 mpfr_t *D, \
			 uint64_t n, uint64_t p, uint64_t q,\
			 mpfr_t *x_low, mpfr_t *x_up,\
			 mpfr_t *y_low, mpfr_t *y_up,\
			 mpfr_t *u_low, mpfr_t *u_up);

void BigMatrixConstruct(mpfr_t *Malade,\
			 mpfr_t *A,mpfr_t *B, mpfr_t *C, mpfr_t *D, mpfr_t *rightSide,\
			 uint64_t size, uint64_t cols, uint64_t n, uint64_t p, uint64_t q);

void LPcomputeRightSide(mpfr_t *rightSide, mpfr_t *A, mpfr_t *B, mpfr_t *C,mpfr_t *D,\
			 uint64_t n, uint64_t p, uint64_t q,\
			 mpfr_t *x_low, mpfr_t *x_up,\
			 mpfr_t *y_low, mpfr_t *y_up,\
			 mpfr_t *u_low, mpfr_t *u_up);

void minOnLine(mpfr_exp_t *v, mpfr_t *A, uint64_t n, uint64_t m);



#endif
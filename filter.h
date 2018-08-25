/*
 ============================================================================
 Name        : filter.h
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : 	  

 ============================================================================
*/

#ifndef _FILTER_H_
#define _FILTER_H_

#include <stdio.h>
#include <stdlib.h>
#include <matrixalg.h>


typedef struct filter_struct
{
	uint64_t n;
	uint64_t p;
	uint64_t q;
	double *A;
	double *B;
	double *C;
	double *D;
	double *ur;		
	uint64_t *w;	
	int *msbU;
	int *lsbU;
}filter;


void filter_allocate			(filter *f, uint64_t n, uint64_t p, uint64_t q);
void filter_free				(filter *f);
void filter_set 				(filter *f, uint64_t n, uint64_t p, uint64_t q, double *A, \
					 								double *B, double *C, double *D, double *ur);
void filter_set_uFormat			(filter *f, int *msb, int *lsb);
void filter_set_w 				(filter *f, uint64_t *w);
void filter_set_w_eql			(uint64_t* wl, uint64_t w_z, uint64_t length);
void filter_set_ur_to_eps		(filter *f, uint64_t *msb, uint64_t *w);
void filter_print				(FILE *file, filter *f);


void createFilterHz1			(filter *Hz1, filter *H);
void createFilterDelta			(filter *Delta, filter *H, uint64_t *msb, uint64_t *w);




#endif

/*
 ============================================================================
 Name        : FXPF.c
 Author      : Anastasia Lozanova Volkova
 Version     : 0.1
 Copyright   : 
 Description : 	This program launches an example of using FXPF library

 ============================================================================
 */

 #include <stdio.h>
 #include "fxpf.h"

int main(int argc, char *argv[] ){


	FILE *input = fopen("./examples/simple.txt", "r");

	int wordlength = 16;
	int wordlength_min = 4;
	if (argc > 2)
	{
		wordlength = atoi(argv[1]);
		wordlength_min = atoi(argv[2]);
		printf("User chose wodlengths between %d and %d \n", wordlength, wordlength_min);
	}

 	// FILE *input = fopen("./examples/input3_f53bits.txt", "r");
	// FILE *input = fopen("./examples/input2.txt", "r");
	// FILE *output = fopen("./examples/output.txt", "a+");
	// FILE *input = fopen("./examples/DL_Cor.txt", "r");
	// FILE *output = fopen("./examples/output_RandSS_new.txt", "a+");
	// FILE *input = fopen("./examples/RandSS.txt", "r");
	// FILE *output = fopen("./examples/output_test2.txt", "a+");
	// FILE *input = fopen("./examples/test2.txt", "r");
	// FILE *output = fopen("./examples/output_test_slides.txt", "a+");
	//FILE *input = fopen("./examples/test_slides.txt", "r");
	//FILE *output = fopen("./examples/output_article.txt", "a+");
	//FILE *input = fopen("./examples/test_article.txt", "r");
	//FILE *input = fopen("./examples/test_thesis6.txt", "r");

	
	//FILE *input = fopen("./examples/thesis_nasty.txt", "r");
	//FILE *input = fopen("./examples/xilinx.txt", "r");

	FILE *output = fopen("./examples/output_RAIM.txt", "a+");
	//FILE *input = fopen("./examples/test_RAIM.txt", "r");
	//FILE *inputU = fopen("./examples/RAIM_input1.txt", "r");

 	if(input == NULL)
 	{
 		fprintf(stderr, "Error opening example file \n");
 		return -1;
 	}

 	if(output == NULL)
 	{
 		fprintf(stderr, "Error opening output for the example file \n");
 		return -1;
 	}

 	uint64_t tmp, tmp2;

	filter H;
	int filterNo = 0;
	double *A;
	double *B;
	double *C;
	double *D;
	double *u_bound;
	uint64_t *wl;
	uint64_t *msb;
	int *lsb;

 	while(!feof(input) && filterNo < 1)
 	{
 	
 		filterNo++;
	 	/* Reading matrices */
		uint64_t n;	
		fscanf(input,"%llu %llu \n", &n, &n);	
		
		A = (double*)malloc(n*n*sizeof(double));
		readDoubleMatrix(input, A,n, n);
		
		uint64_t q;
		fscanf(input,"%llu %llu \n", &tmp, &q);	
		if(tmp != n) 
		{
			printf("Error, incorrect input file format, matrix sizes are not coherent. Must have %llu x %llu, but have %llu x %llu \n", n, q, tmp, q);
			return -1;
		}
		tmp = 0;
		B = (double*)malloc(n*q*sizeof(double));
		readDoubleMatrix(input, B,n, q);
		
		uint64_t p;
		fscanf(input,"%llu %llu \n", &p, &tmp);	
		if(tmp != n) 
		{
			printf("Error, incorrect input file format, matrix sizes are not coherent C. \n");
			return -1;
		}
		tmp = 0;
		C = (double*)malloc(p*n*sizeof(double));
		readDoubleMatrix(input, C,p, n);
		
		
		fscanf(input,"%llu %llu \n", &tmp, &tmp2);
		if(tmp != p || tmp2 != q) 
		{
			printf("Error, incorrect input file format, matrix sizes are not coherent D. \n");
			return -1;
		}
		tmp = 0; tmp2 = 0;	
		D = (double*)malloc(p*q*sizeof(double));
		readDoubleMatrix(input, D,p, q);
		
		fscanf(input,"%llu %llu \n", &tmp, &tmp2);
		if(tmp2 != q || tmp != 1) 
		{
			printf("Error, incorrect input file format, matrix sizes are not coherent u. \n");
			return -1;
		}
		tmp = 0; tmp2 = 0;	

		u_bound = (double*)malloc(1*q*sizeof(double));
		readDoubleMatrix(input, u_bound, 1, q);



/*
		 int *msbU = (int*)malloc(q * sizeof(int));
		 int *lsbU = (int*)malloc(q * sizeof(int));
		 fscanf(inputU,"%llu %llu \n", &tmp, &tmp);
		 readINTMatrix(input, msbU, q, 1);
		 fscanf(inputU,"%llu %llu \n", &tmp, &tmp);
		 readINTMatrix(inputU, lsbU, q, 1);

		 uint64_t T;
		 fscanf(input,"%llu %llu \n", &T, &tmp);
		 double *Ud = (double*)malloc(tmp*T*sizeof(double));
		 readDoubleMatrix(input, Ud, T, tmp);

		
		 mpfr_t *U;
		 U = allocateMPFRMatrix(T, tmp, 64);
		 doubleToMPFRMatrix(U, Ud, T, tmp);
		*/
		
		filter_allocate(&H, n, p, q);
		filter_set(&H, n, p, q, A, B, C, D, u_bound);
		wl = (uint64_t*)calloc(H.p + H.n, (H.p + H.n) * sizeof(uint64_t));
		
		 printf("The filter:\n");
		 filter_print(stderr, &H);

		msb = (uint64_t*)calloc(H.p + H.n, (H.p + H.n) * sizeof(uint64_t));
		lsb = (int*)calloc(H.p + H.n, (H.p + H.n) * sizeof(int));

		

		fxpf_result result;
		fxpf_result_allocate(&result,H.p + H.n);

        /*
		while(wordlength >= wordlength_min)
		{
			filter_set_w_eql(wl, (uint64_t)wordlength, H.p + H.n);

			if(determineFXPF(msb, lsb, &H, wl, &result) < 0)
			{
				fprintf(stderr, "Error determining Fixed-Point Formats.\n");
				break;
			}
			else
			{
				printf("########## FINAL Result #############################################\n");
				printf("Filter #%d: %d additional steps \n", filterNo, result.additionalSteps);
				printf("--------------------------------------------------------------\n");
				printf("Wordlengths: \n");
				UINT64MatrixPrint(stderr, wl, 1, H.p + H.n);
				printf("MSBs: \n");
				UINT64MatrixPrint(stderr, msb, 1, H.p + H.n);
				printf("LSBs: \n");
				INTMatrixPrint(stderr, lsb, 1, H.p + H.n);
				printf("####################################################################\n");

				if(result.additionalSteps > 1)
				{
					fprintf(output, "This filter had %d additional steps with w= %d\n",result.additionalSteps, wordlength);
					filter_print(output, &H);
					fxpf_result_print(output, &result);
					// scanf("lolo");
				}




	
			}
			wordlength--;
			//("%d", &wordlength);
		}
		*/

		filter_set_w_eql(wl, (uint64_t)wordlength, H.p + H.n);
		printf("---------------computing using the DetermineFXPF function...\n");
	    if(determineFXPF(msb, lsb, &H, wl, &result) < 0)
			{
				fprintf(stderr, "Error determining Fixed-Point Formats.\n");
				break;
			}
		else
		{
			printf("########## FINAL Result #############################################\n");
			printf("Filter #%d: %d additional steps \n", filterNo, result.additionalSteps);
			printf("--------------------------------------------------------------\n");
			printf("Wordlengths: \n");
    		UINT64MatrixPrint(stderr, wl, 1, H.p + H.n);
			printf("MSBs: \n");
			UINT64MatrixPrint(stderr, msb, 1, H.p + H.n);
			printf("LSBs: \n");
			INTMatrixPrint(stderr, lsb, 1, H.p + H.n);
			printf("####################################################################\n");

		}

		printf("---------------computing using the FXPF function...\n");
		//if(determineFXPF(msb, lsb, &H, wl, &result) < 0)
//		int FXPF(uint64_t *msb, int *lsb, int *additionalSteps, double *error,
//int *wl, double *A, double *B, double *C, double *D, int n, int p, int q, double *u_bound)

        double *error = (double*)malloc((p + n) * sizeof(double));
        int additionalSteps = -1;
		if(FXPF(msb, lsb, &additionalSteps, error, wl, A, B, C, D, n, p, q, u_bound))
		{
			fprintf(stderr, "Error determining Fixed-Point Formats.\n");
			break;
		}
	    else
	    {
            printf("\n\n\n--------- FXPF Result ----------\n");
		    printf("Additional steps: %d,\n", additionalSteps);
		    printf("--------------------------------------------------------------\n");
		    printf("Wordlengths: \n");
		    UINT64MatrixPrint(stderr, wl, 1, p + n);
		    printf("MSBs: \n");
		    UINT64MatrixPrint(stderr, msb, 1, p + n);
		    printf("LSBs: \n");
		    INTMatrixPrint(stderr, lsb, 1, p + n);
		    printf("Errors: \n");
		    DoubleMatrixPrint(stderr, error, 1, p + n);
		    printf("####################################################################\n");
	    }




		filter_free(&H);
		fxpf_result_free(&result);
		free(u_bound);
		free(wl);
		free(lsb);
		free(msb);
		free(A);
		free(B);
		free(C);
		free(D);

	}
	

	

	
	fclose(input);
	fclose(output);
	return 0;
}





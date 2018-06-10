#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>


void malloc_matrix(int m, int n, double ***matptr);

void multiply(int m1, int n1, double **mat1, int m2, int n2, double **mat2, double **mat3);

void multiply2(int m1, int n1, double **mat1, int m2, int n2, double **mat2, double **mat3);

#define INPUTFILE "input/test2"


int main(int argc,char*argv []) {
	int flag = 0;
	if(argc>1){
		if(atoi(argv[1])==1)flag = 1;
		
		else if(atoi(argv[1])==2)flag=2;
	}
	unsigned ma = 0;	
	unsigned na = 0;
	unsigned mb = 0;	
	unsigned nb = 0;
	int i = 0; // used in loop
	int j = 0; // used in loop
	double **A, **B, **C ,**C2;
	FILE *infile, *outfile;

	
	infile = fopen(INPUTFILE, "r");
	if(infile == NULL) {
    	printf("Error in Opening infile");
    	return EXIT_FAILURE;
	}

	fscanf(infile, "%u %u", &ma, &na);	// read A's dimensions
	
	malloc_matrix(ma, na, &A);
	
	for(i = 0; i < ma; i++)	// read A's content
		for(j = 0; j < na; j++)
			fscanf(infile, "%lf", &A[i][j]);
	
	fscanf(infile, "%u %u", &mb, &nb);	// read B's dimensions
	
	if(na != mb) {	// check dimension
		printf("matrices dimension error \n");		
		exit(1);	
	}

	malloc_matrix(mb, nb, &B);
	
	for(i = 0; i < mb; i++)	// read B's content
		for(j = 0; j < nb; j++)
			fscanf(infile, "%lf", &B[i][j]);
	fclose(infile);
	
	
	malloc_matrix(ma, nb, &C);
	malloc_matrix(ma, nb, &C2);
	/*** do multiplication in different methods ***/
	outfile = fopen("output_orig.txt", "w+");
	
	double st,en;
	if(flag!=2){
	st=omp_get_wtime();
	multiply2(ma, na, A, mb, nb, B, C2);
	en=omp_get_wtime();
	printf("Column major: %lf\n",en-st);	
	fprintf(outfile, "Column major: %lf s \n", en-st);
	}
	
	// method 2 : row major
	if(flag!=1){
	st=omp_get_wtime();
	multiply(ma, na, A, mb, nb, B, C);
	en=omp_get_wtime();
	printf("Row major: %lf\n",en-st);	

	fprintf(outfile, "Row major: %lf s \n", en-st);
	}
	// method 2 end
	/*
	for(i = 0; i < ma; i++) 
		for(j = 0; j < nb; j++)
			if(C2[i][j]!=C[i][j])printf("error\n");
	*/
	// output the result matrix C 
	fprintf(outfile, "%d %d \n", ma, nb);
	
	for(i = 0; i < ma; i++) {
		for(j = 0; j < nb; j++)
			if(flag!=1)
				fprintf(outfile, "%8.0lf ", C[i][j]);
			else if(flag!=2)
				fprintf(outfile, "%8.0lf ", C2[i][j]);
		fprintf(outfile, "\n");	
	}
	fclose(outfile);

	return 0;
}

void malloc_matrix(int m, int n, double ***matptr) {
	/*
	// the memory of the 2d array is not consecutive (each row)
	int i;	
	*matptr = malloc(m * sizeof(double *));
	for(i = 0; i < m; i++) 
		(*matptr)[i] = malloc(n * sizeof(float));
	*/
	
	// the memory of the 2d array is consecutive (each row)
	int i;
	double *tmp;	
	tmp = malloc(m * n * sizeof(double *));
	*matptr = malloc(m * sizeof(double *));
	for(i = 0; i < m; i++)
		(*matptr)[i] = &(tmp[n*i]);
	
}

void multiply(int m1, int n1, double **mat1, int m2, int n2, double **mat2, double **mat3) {	
	int i, j, k;
	for(i = 0; i < m1; i++)
		for(j = 0; j < n2; j++)
			mat3[i][j] = 0;	
	for(i = 0; i < m1; i++)
		for(j = 0; j < n1; j++)
			for(k = 0; k < n2; k++)
				mat3[i][k] += mat1[i][j] * mat2[j][k];
}

void multiply2(int m1, int n1, double **mat1, int m2, int n2, double **mat2, double **mat3) {	
	int i, j, k;
	for(i = 0; i < m1; i++)
		for(j = 0; j < n2; j++)
			mat3[i][j] = 0;	
	for(i = 0; i < m1; i++)
		for(j = 0; j < n2; j++)
			for(k = 0; k < n1; k++)
				mat3[i][j] += mat1[i][k] * mat2[k][j];
}



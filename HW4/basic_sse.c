#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <x86intrin.h>
#include <immintrin.h>

#define INPUTFILE "input/test2"

void malloc_matrix(int m, int n, double ***matptr);
void multiply(int m1, int n1, double **mat1, int m2, int n2, double **mat2, double **mat3);
void multiply_NotP(int m1, int n1, double **mat1, int m2, int n2, double **mat2,double **mat3);

int to_x_multiple(int x, int in);
int thread_count=1;
double **A, **B, **C;

int main(int argc,char*argv []) {
	//輸入的
	unsigned ma_in = 0;	
	unsigned na_in = 0;
	unsigned mb_in = 0;	
	unsigned nb_in = 0;
	//真的allocate的size
	unsigned ma = 0;	
	unsigned na = 0;
	unsigned mb = 0;	
	unsigned nb = 0;

	int i = 0; // used in loop
	int j = 0; // used in loop

	FILE *infile, *outfile;
	
	infile = fopen(INPUTFILE, "r");
	if(infile == NULL) {
    	printf("Error in Opening infile");
    	return EXIT_FAILURE;
	}

	fscanf(infile, "%u %u", &ma_in, &na_in);	// read A's dimensions
	//補滿4
	ma = to_x_multiple(4, ma_in);
	na = to_x_multiple(4, na_in);
	
	malloc_matrix(ma, na, &A);
	
	for(i = 0; i < ma_in; i++)	{
		for(j = 0; j < na_in; j++) {
				fscanf(infile, "%lf", &A[i][j]);
		}
	}

	fscanf(infile, "%u %u", &mb_in, &nb_in);	// read B's dimensions
	mb = to_x_multiple(4, mb_in);
	nb = to_x_multiple(4, nb_in);
	
	if(na != mb) {	// check dimension
		printf("matrices dimension error %d %d %d %d \n",na,nb,ma,mb);		
		exit(1);	
	}

	
     	malloc_matrix(mb, nb, &B);
	for(i = 0; i < mb_in; i++)	{
		for(j = 0; j < nb_in; j++) {
				fscanf(infile, "%lf", &B[i][j]);
		}
	}
	
	fclose(infile);
	
	malloc_matrix(ma, nb, &C);
	
	/*** do multiplication in different methods ***/
	outfile = fopen("output_sseOpenmp.txt", "w+");
	
	// method 1 : traditional
	thread_count = strtol(argv[1],NULL,10);

	//not parallel

	double st=omp_get_wtime();
	multiply_NotP(ma, na, A, mb, nb, B, C);
	double en=omp_get_wtime();
	printf("Basic + SSE: %lf\n",en-st);	
	fprintf(outfile, "Basic + SSE: %lf s \n", en-st);

	//parallel

	st=omp_get_wtime();
	multiply(ma, na, A, mb, nb, B, C);
	en=omp_get_wtime();
	printf("Basic + SSE + openmp: %lf\n",en-st);

	fprintf(outfile, "Basic + SSE + openmp: %lf s \n", en-st);
	// method 1 end

	// output the result matrix C 
	fprintf(outfile, "%d %d \n", ma, nb);
	for(i = 0; i < ma; i++) {
		for(j = 0; j < nb; j++)
			fprintf(outfile, "%8.0lf ", C[i][j]);
		fprintf(outfile, "\n");	
	}
	fclose(outfile);

	return 0;
}

int to_x_multiple(int x, int in) {
	if(in%x != 0)
		in += 4 - in%4;
	return in;
}

void malloc_matrix(int m, int n, double ***matptr) {
	// the memory of the 2d array is consecutive (each row)
	int i;
	double *tmp;	
	tmp = _mm_malloc(m * n * sizeof(double), 64);
//	tmp = malloc(m * n * sizeof(double));
	memset(tmp, 0, m*n * sizeof (double));
	*matptr = malloc(m * sizeof(double *));
	for(i = 0; i < m; i++)
		(*matptr)[i] = &(tmp[n*i]);
	
}


void multiply(int m1, int n1, double **mat1, int m2, int n2, double **mat2, double **mat3) {	
	int i, j, k;

#pragma omp parallel for num_threads(thread_count)private(i,j) shared(mat3)
	for(i = 0; i < m1; i++)
		for(j = 0; j < n2; j++)
			mat3[i][j] = 0;	
	#pragma omp parallel for num_threads(thread_count) private(i, j, k) shared(mat3,mat1,mat2,m1,n1,n2)
    for(i=0; i<m1; i++) {
        for(j=0; j<n1; j++) {
            __m256d a4 = _mm256_set_pd(mat1[i][j],mat1[i][j],mat1[i][j],mat1[i][j]);
            for(k=0; k<n2; k+=4) {
                __m256d c4 = _mm256_load_pd(&mat3[i][k]);
                __m256d b4 = _mm256_load_pd(&mat2[j][k]);
                c4 = _mm256_add_pd(_mm256_mul_pd(a4,b4),c4);
                _mm256_store_pd(&mat3[i][k], c4);
            }
        }
    }
}

void multiply_NotP(int m1, int n1, double **mat1, int m2, int n2, double **mat2, double **mat3) {	
	int i, j, k;

	for(i = 0; i < m1; i++)
		for(j = 0; j < n2; j++)
			mat3[i][j] = 0;	
	
	for(i=0; i<m1; i++) {
        for(j=0; j<n1; j++) {
            __m256d a4 = _mm256_set_pd(mat1[i][j],mat1[i][j],mat1[i][j],mat1[i][j]);
            for(k=0; k<n2; k+=4) {
                __m256d c4 = _mm256_load_pd(&mat3[i][k]);
                __m256d b4 = _mm256_load_pd(&mat2[j][k]);
                c4 = _mm256_add_pd(_mm256_mul_pd(a4,b4),c4);
                _mm256_store_pd(&mat3[i][k], c4);
            }
        }
    }
}

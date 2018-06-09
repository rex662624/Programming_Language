#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <x86intrin.h>

#define INPUTFILE "input/input1000.txt"

void malloc_matrix(int m, int n, float ***matptr);
void add(int m, int n, float **mat1, float **mat2, float **mat3);
void sub(int m, int n, float **mat1, float **mat2, float **mat3);

// split one matrix to four same-shape sub-matrices
void matrix_split(int m, int n, float **mat,
					float **mat1, float **mat2,
					float **mat3, float **mat4);

// merge four same-shape sub-matrices to one matrix
void matrix_merge(int m, int n, float **mat,
					float **mat1, float **mat2,
					float **mat3, float **mat4);

void multiply(int m1, int n1, float **mat1, int m2, int n2, float **mat2, float **mat3);
void multiply_NotP(int m1, int n1, float **mat1, int m2, int n2, float **mat2,float **mat3);

int to_x_multiple(int x, int in);
int thread_count=1;
float **A, **B, **C;

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

	float **A11, **A12, **A21, **A22;
	float **B11, **B12, **B21, **B22;
	float **C11, **C12, **C21, **C22;
	float **C111, **C112, **C221, **C222;
	float **P1, **P2, **P, **Q1, **Q, **R1, **R, **S1, **S, **T1, **T, **U1, **U2, **U, **V1, **V2, **V;
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
				fscanf(infile, "%f", &A[i][j]);
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
				fscanf(infile, "%f", &B[i][j]);
		}
	}
	
	fclose(infile);
	
	malloc_matrix(ma, nb, &C);
	
	/*** do multiplication in different methods ***/
	outfile = fopen("output_sseOpenmp.txt", "w+");
/*	
	// method 1 : traditional
	thread_count = strtol(argv[1],NULL,10);
	double st,en;
	//not parallel
	double st=omp_get_wtime();
	multiply_NotP(ma, na, A, mb, nb, B, C);
	double en=omp_get_wtime();
	printf("Basic: %lf\n",en-st);	
	fprintf(outfile, "Basic: %lf s \n", en-st);
	//parallel
	st=omp_get_wtime();
	multiply(ma, na, A, mb, nb, B, C);
	en=omp_get_wtime();
	printf("Basic+SSE +openmp: %lf\n",en-st);

	fprintf(outfile, "Basic+SSE +openmp: %lf s \n", en-st);
	// method 1 end
*/
	double st,en;
	// method 2 : Strassen algorithm
	st=omp_get_wtime();
	int na_2,nb_2,ma_2,mb_2;
	//需要做偶數check
	na_2= to_x_multiple(4, na/2);
	nb_2= to_x_multiple(4, nb/2);
	ma_2= to_x_multiple(4, ma/2);
	mb_2= to_x_multiple(4, nb/2);
	// malloc
	malloc_matrix(ma_2, na_2, &A11);
	malloc_matrix(ma_2, na_2, &A12);
	malloc_matrix(ma_2, na_2, &A21);
	malloc_matrix(ma_2, na_2, &A22);
	malloc_matrix(mb_2, nb_2, &B11);
	malloc_matrix(mb_2, nb_2, &B12);
	malloc_matrix(mb_2, nb_2, &B21);
	malloc_matrix(mb_2, nb_2, &B22);
	malloc_matrix(ma_2, nb_2, &C11);
	malloc_matrix(ma_2, nb_2, &C12);
	malloc_matrix(ma_2, nb_2, &C21);
	malloc_matrix(ma_2, nb_2, &C22);
	
	malloc_matrix(ma_2, nb_2, &C111);
	malloc_matrix(ma_2, nb_2, &C112);
	malloc_matrix(ma_2, nb_2, &C221);
	malloc_matrix(ma_2, nb_2, &C222);

	malloc_matrix(ma_2, na_2, &P1);
	malloc_matrix(mb_2, nb_2, &P2);
	malloc_matrix(ma_2, nb_2, &P);

	malloc_matrix(ma_2, na_2, &Q1);
	malloc_matrix(ma_2, nb_2, &Q);

	malloc_matrix(mb_2, nb_2, &R1);
	malloc_matrix(ma_2, nb_2, &R);

	malloc_matrix(mb_2, nb_2, &S1);
	malloc_matrix(ma_2, nb_2, &S);

	malloc_matrix(ma_2, na_2, &T1);
	malloc_matrix(ma_2, nb_2, &T);

	malloc_matrix(ma_2, na_2, &U1);
	malloc_matrix(mb_2, nb_2, &U2);
	malloc_matrix(ma_2, nb_2, &U);

	malloc_matrix(ma_2, na_2, &V1);
	malloc_matrix(mb_2, nb_2, &V2);
	malloc_matrix(ma_2, nb_2, &V);

	// calculate (may do parallel below)
	
	matrix_split(ma, na, A, A11, A12, A21, A22);
	matrix_split(mb, nb, B, B11, B12, B21, B22);
	
	thread_count = strtol(argv[1],NULL,10);

		add(ma/2, na/2, A11, A22, P1); // P
		add(mb/2, nb/2, B11, B22, P2);

		multiply(ma/2, na/2, P1, mb/2, nb/2, P2, P);
		
		add(ma/2, na/2, A21, A22, Q1); // Q
		multiply(ma/2, na/2, Q1, mb/2, nb/2, B11, Q);
		
		sub(mb/2, nb/2, B12, B22, R1); // R
		multiply(ma/2, na/2, A11, mb/2, nb/2, R1, R);
		


		sub(mb/2, nb/2, B21, B11, S1); // S
		multiply(ma/2, na/2, A22, mb/2, nb/2, S1, S);


		add(ma/2, na/2, A11, A12, T1); // T
		multiply(ma/2, na/2, T1, mb/2, nb/2, B22, T);
	
		sub(ma/2, na/2, A21, A11, U1); // U
		add(mb/2, nb/2, B11, B12, U2);
		multiply(ma/2, na/2, U1, mb/2, nb/2, U2, U);
	
		sub(ma/2, na/2, A12, A22, V1); // V
		add(mb/2, nb/2, B21, B22, V2);
		
		multiply(ma/2, na/2, V1, mb/2, nb/2,V2,V);
	
		add(ma/2, nb/2, P, S, C111); // C11
		sub(ma/2, nb/2, C111, T, C112);
		add(ma/2, nb/2, C112, V, C11);
		
		add(ma/2, nb/2, R, T, C12); // C12
	
		add(ma/2, nb/2, Q, S, C21); // C21
	
		add(ma/2, nb/2, P, R, C221); // C22
		sub(ma/2, nb/2, C221, Q, C222);
		add(ma/2, nb/2, C222, U, C22);
	

	// merge to C
	matrix_merge(ma, nb, C, C11, C12, C21, C22);

	en=omp_get_wtime();
 	printf("Strassen 1 cut + openpmp + SSE: %lf\n",en-st);
	fprintf(outfile, "Strassen 1 cut + openpmp + SSE %lf s \n", en-st);
	// method 2 end
	

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

void malloc_matrix(int m, int n, float ***matptr) {
	// the memory of the 2d array is consecutive (each row)
	int i;
	float *tmp;	
	tmp = _mm_malloc(m * n * sizeof(float), 16);
//	tmp = malloc(m * n * sizeof(float));
	memset(tmp, 0, m*n * sizeof (float));
	*matptr = malloc(m * sizeof(float *));
	for(i = 0; i < m; i++)
		(*matptr)[i] = &(tmp[n*i]);
	
}

void add(int m, int n, float **mat1, float **mat2, float **mat3) {	
	int i, j;
	#pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat3,mat2,mat1)
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			mat3[i][j] = mat1[i][j] + mat2[i][j];	
}

void sub(int m, int n, float **mat1, float **mat2, float **mat3) {	
	int i, j;
	#pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat3,mat2,mat1)
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			mat3[i][j] = mat1[i][j] - mat2[i][j];	
}

void multiply(int m1, int n1, float **mat1, int m2, int n2, float **mat2, float **mat3) {	
	int i, j, k;

#pragma omp parallel for num_threads(thread_count)private(i,j) shared(mat3)
	for(i = 0; i < m1; i++)
		for(j = 0; j < n2; j++)
			mat3[i][j] = 0;	
	#pragma omp parallel for num_threads(thread_count) private(i, j, k) shared(mat3,mat1,mat2,m1,n1,n2)
    for(i=0; i<m1; i++) {
        for(j=0; j<n1; j++) {
            __m128 a4 = _mm_load1_ps(&mat1[i][j]);
            for(k=0; k<n2; k+=4) {
                __m128 c4 = _mm_load_ps(&mat3[i][k]);
                __m128 b4 = _mm_load_ps(&mat2[j][k]);
                c4 = _mm_add_ps(_mm_mul_ps(a4,b4),c4);
                _mm_store_ps(&mat3[i][k], c4);
            }
        }
    }
}
void multiply_NotP(int m1, int n1, float **mat1, int m2, int n2, float **mat2, float **mat3) {	
	int i, j, k;

	for(i = 0; i < m1; i++)
		for(j = 0; j < n2; j++)
			mat3[i][j] = 0;	

	for(i = 0; i < m1; i++)
		for(j = 0; j < n1; j++)
			for(k = 0; k < n2; k++)
				mat3[i][k] += mat1[i][j] * mat2[j][k];
}				

void matrix_split(int m, int n, float **mat,
					float **mat11, float **mat12,
					float **mat21, float **mat22) {
	int i, j;
	if(m % 2 != 0 || n % 2 != 0) {
		printf("m = %d, n = %d \n", m, n);
		printf("error, cannot be split to 4 sub-matrices \n");
		exit(1);
	}
	
#pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat,mat11,m,n)	
	for(i = 0; i < m/2; i++)
		for(j = 0; j < n/2; j++)
			mat11[i][j] = mat[i][j];

 #pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat,mat12,m,n)
	for(i = 0; i < m/2; i++)
		for(j = n/2; j < n; j++)
			mat12[i][j - n/2] = mat[i][j];
 #pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat,mat21,m,n)
	for(i = m/2; i < m; i++)
		for(j = 0; j < n/2; j++)
			mat21[i - m/2][j] = mat[i][j];

 #pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat,mat22,m,n)
	for(i = m/2; i < m; i++)
		for(j = n/2; j < n; j++)
			mat22[i - m/2][j - n/2] = mat[i][j];
	
}

void matrix_merge(int m, int n, float **mat,
					float **mat11, float **mat12,
					float **mat21, float **mat22) {
	int i, j;
	if(m % 2 != 0 || n % 2 != 0) {
		printf("error, cannot merge by 4 sub-matrices \n");
		exit(1);
	}
 #pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat,mat11,n,m)
	for(i = 0; i < m/2; i++)
		for(j = 0; j < n/2; j++)
			mat[i][j] = mat11[i][j];

 #pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat,mat12,n,m)
	for(i = 0; i < m/2; i++)
		for(j = n/2; j < n; j++)
			mat[i][j] = mat12[i][j - n/2];
#pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat,mat21,n,m)      
	for(i = m/2; i < m; i++)
		for(j = 0; j < n/2; j++)
			mat[i][j] = mat21[i - m/2][j];
#pragma omp parallel for num_threads(thread_count) private(i,j) shared(mat,mat22,n,m)      
	for(i = m/2; i < m; i++)
		for(j = n/2; j < n; j++)
			mat[i][j] = mat22[i - m/2][j - n/2];
}

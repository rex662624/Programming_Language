//c99 -O3 -fopenmp -Wall foo.c
#include <stdio.h>
#include <string.h>
#include <x86intrin.h>
#include <omp.h>

#define INPUTFILE "input1000.txt"

void malloc_matrix(int m, int n, float ***matptr) {
	// the memory of the 2d array is consecutive (each row)
	int i;
	float *tmp;	
//	tmp = _mm_malloc(m * n * sizeof(float), 16);
	tmp = malloc(m * n * sizeof(float));
	memset(tmp, 0, m*n * sizeof (float));
	*matptr = malloc(m * sizeof(float *));
	for(i = 0; i < m; i++)
		(*matptr)[i] = &(tmp[n*i]);
	
}

void gemm(float ** restrict a, float ** restrict b, float ** restrict c, int ma, int na, int nb) {
	int i, j, k;
    for(i = 0; i < ma; i++) {
        for(j = 0; j < na; j++) {
            for(k = 0; k < nb; k++) {
                c[i][k] += a[i][ j] * b[j][k];
            }
        }
    }
}

void gemm_tlp(float ** restrict a, float ** restrict b, float ** restrict c, int ma, int na, int nb) {
	int i, j, k;
	#pragma omp parallel for private(i, j, k) shared(a, b, c)
    for(i=0; i<ma; i++) {
        for(j=0; j<na; j++) {
            for(k=0; k<nb; k++) {
				c[i][k] += a[i][j] * b[j][k];
            }
        }
    }
}   

void gemm_tlp_simd(float ** restrict a, float ** restrict b, float ** restrict c, int ma, int na, int nb) {
	int i, j, k;
	#pragma omp parallel for private(i, j, k) shared(a, b, c)
    for(i=0; i<ma; i++) {
        for(j=0; j<na; j++) {
            __m128 a4 = _mm_set1_ps(a[i][j]);
            for(k=0; k<nb; k+=4) {
                __m128 c4 = _mm_load_ps(&c[i][k]);
                __m128 b4 = _mm_load_ps(&b[j][k]);
                c4 = _mm_add_ps(_mm_mul_ps(a4,b4),c4);
                _mm_store_ps(&c[i][k], c4);
            }
        }
    }
}

int to_x_multiple(int x, int in) {
	if(in%x != 0)
		in += 4 - in%4;
	return in;
}


int main(void) {
	//輸入的
	int ma_in = 0;	
	int na_in = 0;
	int mb_in = 0;	
	int nb_in = 0;
	//真的allocate的size
	int ma = 0;
    int na = 0;
	int mb = 0;
	int nb = 0;
	int i, j, k;

	FILE *infile, *outfile;
	infile = fopen(INPUTFILE, "r");
	if(infile == NULL) {
    	printf("Error in opening infile");
    	return EXIT_FAILURE;
	}

	fscanf(infile, "%d %d", &ma_in, &na_in);
	ma = to_x_multiple(4, ma_in);
	na = to_x_multiple(4, na_in);
    //declare array a  (float *a = _mm_malloc(ma*na * sizeof *a, 16);)
	float **a;
     	malloc_matrix(ma, na, &a);
	for(i = 0; i < ma; i++)	{
		for(j = 0; j < na; j++) {
			if(i < ma_in && j < na_in)
				fscanf(infile, "%f", &a[i][j]);
			// this line became redundant becasue malloc_matrix already initilize
			/*else
				a[i][j] = 0;*/
		}
	}
	fscanf(infile, "%d %d", &mb_in, &nb_in);
	mb = to_x_multiple(4, mb_in);
	nb = to_x_multiple(4, nb_in);
	//declare array b (float *b = _mm_malloc(na*nb * sizeof *b, 16);)
	float **b;
     	malloc_matrix(mb, nb, &b);
	for(i = 0; i < mb; i++)	{
		for(j = 0; j < nb; j++) {
			if(i < mb_in && j < nb_in)
				fscanf(infile, "%f", &b[i][j]);
			else
				b[i][j] = 0;
		}
	}
	fclose(infile);
    /*
    float *c1 = _mm_malloc(ma*nb * sizeof *c1, 16);
    float *c2 = _mm_malloc(ma*nb * sizeof *c2, 16);
    float *c3 = _mm_malloc(ma*nb * sizeof *c2, 16);
    memset(c1, 0, ma*nb * sizeof *c1);
    memset(c2, 0, ma*nb * sizeof *c2);
    memset(c3, 0, ma*nb * sizeof *c3);
     */
    float **c1,**c2,**c3;
   malloc_matrix(ma, nb, &c1);	
   malloc_matrix(ma, nb, &c2);	
   malloc_matrix(ma, nb, &c3);	  
  
   double dtime;
    dtime = -omp_get_wtime();
    gemm(a,b,c1,ma, na, nb);
    dtime += omp_get_wtime();
    printf("time %f\n", dtime);

    dtime = -omp_get_wtime();
    gemm_tlp(a,b,c2,ma, na, nb);
    dtime += omp_get_wtime();
    printf("time %f\n", dtime);

    dtime = -omp_get_wtime();
    gemm_tlp_simd(a,b,c3,ma, na, nb);
    dtime += omp_get_wtime();
     printf("time %f\n", dtime);
//----------------------------------------------------------printf result

	for(i = 0; i < ma; i++) {
		for(j = 0; j < nb; j++) {
			if(i < ma_in && j < nb_in)			
				if(c1[i][j]!=c2[i][j]||c3[i][j]!=c1[i][j])printf("error\n");//printf("%.2f  ", c1[i][j]);
		}		

	}	
/*	printf("\n");

	for(i = 0; i < ma; i++) {
		for(j = 0; j < nb; j++) {
			if(i < ma_in && j < nb_in)			
				printf("%.2f  ", c2[i][ j]);
		}		
		printf("\n");	
	}	
	printf("\n");
*/
/*	for(i = 0; i < ma; i++) {
		for(j = 0; j < nb; j++) {
			if(i < ma_in && j < nb_in)			
				printf("%.2f  ", c3[i*nb + j]);
		}		
		printf("\n");	
	}	
	printf("\n");
*/
	printf("%.2f  ", c1[0][6]);printf("%.2f  ", c2[1][5]);printf("%.2f  ", c3[1][8]);printf("\n");	
    //printf("error %d\n", memcmp(&c1[0][0],&c2[0][0], ma*nb*sizeof *c1));
   //printf("error %d\n", memcmp(&c1[0][0],&c3[0][0], ma*nb*sizeof *c1));
}

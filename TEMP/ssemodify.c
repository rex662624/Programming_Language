//c99 -O3 -fopenmp -Wall foo.c
#include <stdio.h>
#include <string.h>
#include <x86intrin.h>
#include <omp.h>

void gemm(float * restrict a, float * restrict b, float * restrict c, int n) {
    for(int i=0; i<n; i++) {
        for(int k=0; k<n; k++) {
            for(int j=0; j<n; j++) {
                c[i*n+j] += a[i*n+k]*b[k*n+j];
            }
        }
    }
}

void gemm_tlp(float * restrict a, float * restrict b, float * restrict c, int n) {
    #pragma omp parallel for
    for(int i=0; i<n; i++) {
        for(int k=0; k<n; k++) {
            for(int j=0; j<n; j++) {
                c[i*n+j] += a[i*n+k]*b[k*n+j];
            }
        }
    }
}   

void gemm_tlp_simd(float * restrict a, float * restrict b, float * restrict c, int n) {
    #pragma omp parallel for
    for(int i=0; i<n; i++) {
        for(int k=0; k<n; k++) {
            __m128 a4 = _mm_set1_ps(a[i*n+k]);
            for(int j=0; j<n; j+=4) {//一次load 4 
                __m128 c4 = _mm_load_ps(&c[i*n+j]);
                __m128 b4 = _mm_load_ps(&b[k*n+j]);
                c4 = _mm_add_ps(_mm_mul_ps(a4,b4),c4);
                _mm_store_ps(&c[i*n+j], c4);
            }
        }
    }
}


static inline __m256  gemm_tlp_avx(float * restrict a, float * restrict b, float * restrict c, int n) {
   __m256 result;
    #pragma omp parallel for
    for(int i=0; i<n; i++) {
        for(int k=0; k<n; k++) {
            __m256 a4 = _mm256_set1_ps(a[i*n+k]);
            for(int j=0; j<n; j+=8) {//一次load 4 
                __m256 c4 = _mm256_load_ps(&c[i*n+j]);
                __m256 b4 = _mm256_load_ps(&b[k*n+j]);
                c4 = _mm256_add_ps(_mm256_mul_ps(a4,b4),c4);
                _mm256_store_ps(&c[i*n+j], c4);
            }
        }
    }
    return result;
}

int main(void) {
    int n = 1000;
    float *a = _mm_malloc(n*n * sizeof *a, 16);
    float *b = _mm_malloc(n*n * sizeof *b, 16);
    float *c1 = _mm_malloc(n*n * sizeof *c1, 16);
    float *c2 = _mm_malloc(n*n * sizeof *c2, 16);
    float *c3 = _mm_malloc(n*n * sizeof *c3, 16);
	float *c4 = _mm_malloc(n*n * sizeof *c4, 16);	
	
    for(int i=0; i<n*n; i++) scanf("%f",&a[i]);//a[i] = 1.0*i;
    for(int i=0; i<n*n; i++) scanf("%f",&b[i]);//b[i] = 1.0*i;
    memset(c1, 0, n*n * sizeof *c1);
    memset(c2, 0, n*n * sizeof *c2);
    memset(c3, 0, n*n * sizeof *c3);
    double dtime;

    dtime = -omp_get_wtime();
    gemm(a,b,c1,n);
    dtime += omp_get_wtime();
    printf("time %f\n", dtime);

    dtime = -omp_get_wtime();
    gemm_tlp(a,b,c2,n);
    dtime += omp_get_wtime();
    printf("time %f\n", dtime);

    dtime = -omp_get_wtime();
    gemm_tlp_simd(a,b,c3,n);
    dtime += omp_get_wtime();
    printf("time %f\n", dtime);
	
    dtime = -omp_get_wtime();
    gemm_tlp_avx(a,b,c4,n);
    dtime += omp_get_wtime();
	
    printf("time %f\n", dtime);
    printf("error %d\n", memcmp(c1,c2, n*n*sizeof *c1));
    printf("error %d\n", memcmp(c1,c3, n*n*sizeof *c1));
    printf("error %d\n", memcmp(c1,c4, n*n*sizeof *c1));
    printf("%f %f\n",c3[5],c3[0]);
    printf("%f %f\n",c4[5],c4[0]);
}

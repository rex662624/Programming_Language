#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define MAXTHREADS 4
#define INPUTFILE "input/test2"
int arow,acol,brow,bcol;
double ** a;
double ** b;
double ** c;
void *multiply(void *arg);

int main() {
	double start,end;
	int i,j;
	double time;
	FILE *infile, *outfile;
	infile = fopen(INPUTFILE, "r");
	
	if(infile == NULL) {
    	printf("Error in Opening infile");
    	return EXIT_FAILURE;
	}

	fscanf(infile, "%d %d", &arow, &acol);	
	
	a=(double **) calloc(arow, sizeof(double *));
	for(i=0; i<arow; i++){
		a[i]=(double *) calloc(acol, sizeof(double));
	}
	
	for(i=0;i<arow;i++){
		for(j=0;j<acol;j++){
			fscanf(infile, "%lf", &a[i][j]);	
		}
			
	}
	
	fscanf(infile, "%d %d", &brow, &bcol);
	
	b=(double **)calloc(brow, sizeof(double *));
	
	for(i=0; i<brow; i++){
		b[i]=(double *)calloc(bcol, sizeof(double));
	}
	
	for(i=0;i<brow;i++){
		for(j=0;j<bcol;j++){
			fscanf(infile, "%lf", &b[i][j]);
		}
				
	}
		
	fclose(infile);
	
	c=(double **)calloc(arow, sizeof(double *));
	for(i=0; i<arow; i++){
		c[i]=(double *)calloc(bcol, sizeof(double));
	}
	for(i=0;i<arow;i++){
		for(j=0;j<bcol;j++){
			c[i][j] = 0.0;
		}
	}
	
	outfile = fopen("output_basicpthread.txt", "a+");
	pthread_t pt[MAXTHREADS];
	start = omp_get_wtime();
	for(i=0;i<MAXTHREADS;i++){
		pthread_create(&pt[i],NULL,multiply,(void *)i);			
	}
	for(i=0;i<MAXTHREADS;i++){
		pthread_join(pt[i],NULL);
	}
	end = omp_get_wtime();
	time = ((double) (end - start));
	printf("Basic + pthread: %lf s\n", time);
	fprintf(outfile, "Basic + pthread: %lf\n", time);
	fprintf(outfile, "%d %d \n", arow, bcol);
	for(i = 0; i < arow; i++) {
		for(j = 0; j < bcol; j++)
				fprintf(outfile, "%8.0lf ", c[i][j]);
		
		fprintf(outfile, "\n");		
	}
	fclose(outfile);

	return 0;
}

void *multiply(void *arg){
	int i,j,k,l;
	l = (int)arg;
	int part,start,end;
	part = arow/MAXTHREADS;
	start = l*part;
	end = start + part;
	for(i=start;i<end;i++){
		for(j=0;j<acol;j++){
			for(k=0;k<bcol;k++){
				c[i][k] += a[i][j]*b[j][k];
			}
		}
	}
	return NULL;
}

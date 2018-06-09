#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define MAXTHREADS 4
#define INPUTFILE "input/input4.txt"
int arow,acol,brow,bcol;
float ** a;
float ** b;
float ** c;
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
	
	a=(float **) calloc(arow, sizeof(float *));
	for(i=0; i<arow; i++){
		a[i]=(float *) calloc(acol, sizeof(float));
	}
	
	for(i=0;i<arow;i++){
		for(j=0;j<acol;j++){
			fscanf(infile, "%f", &a[i][j]);	
		}
			
	}
	
	fscanf(infile, "%d %d", &brow, &bcol);
	
	b=(float **)calloc(brow, sizeof(float *));
	
	for(i=0; i<brow; i++){
		b[i]=(float *)calloc(bcol, sizeof(float));
	}
	
	for(i=0;i<brow;i++){
		for(j=0;j<bcol;j++){
			fscanf(infile, "%f", &b[i][j]);
		}
				
	}
		
	fclose(infile);
	
	c=(float **)calloc(arow, sizeof(float *));
	for(i=0; i<arow; i++){
		c[i]=(float *)calloc(bcol, sizeof(float));
	}
	for(i=0;i<arow;i++){
		for(j=0;j<bcol;j++){
			c[i][j] = 0.0;
		}
	}
	
	outfile = fopen("outputbasic.txt", "w+");
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
	printf("Time:%lf ms\n", time);
	fprintf(outfile, "Time:%lf s\n", time);
	fprintf(outfile, "%d %d \n", arow, bcol);
	for(i = 0; i < arow; i++) {
		for(j = 0; j < bcol; j++){
			if(j==0)
				fprintf(outfile, "%lf", c[i][j]);
			else
				fprintf(outfile, " %lf", c[i][j]);
		}
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

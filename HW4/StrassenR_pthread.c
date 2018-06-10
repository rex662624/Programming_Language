/*
a  |  b     e  |  f
c  |  d     g  |  h
p1=a(f-h)
p2=(a+b)h
p3=(c+d)e
p4=d(g-e)
p5=(a+d)(e+h)
p6=(b-d)(g+h)
p7=(a-c)(e+f)
c11=p4+p5-p2+p6
c12=p1+p2
c21=p3+p4
c22=p1+p5-p3-p7
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<pthread.h>
#include<omp.h>
#define INPUTFILE "input/input256.txt"

//Structure to hold two matrices and their dimension
typedef struct{
	int ** a;
	int ** b;
	int arow,acol,brow,bcol,size; 
} abmatrix;

//Structure to hold a single matrix and its dimension
typedef struct{
	int ** m;
	int row,col,size;
} onematrix;


int squarematrix(int num){
	int original_num = num, lower_power = 0, i, actual_num = 1;
	if(num==1){ 
		return 1;
	} 
	while(num>1){ 
		lower_power++;
		num /= 2;
	}
	for(i=0;i<lower_power;i++) {
		actual_num *= 2;
	}

	if(actual_num == original_num){
		return original_num;
	} 
	else{
		return actual_num * 2;
	} 
}

int ** creatematrix(int row,int col){
	int i;
	int ** m = (int **) calloc(row, sizeof(int *));
	for(i=0; i<row; i++)
	{
		m[i] = (int *) calloc(col, sizeof(int));
	}
	return m;
}


//function to allocate memory to the two-matrix structure defined above
abmatrix * initstruct2(int resize){
	abmatrix * m = (abmatrix *) malloc(sizeof(abmatrix));
	m->a = creatematrix(resize,resize);
	m->b = creatematrix(resize,resize);
	m->size = resize;
	return m;
}

//function to allocate memory to the one-matrix structure defined above
onematrix * initstruct(int resize){
	onematrix * M = (onematrix *) malloc(sizeof(onematrix));
	M->m = creatematrix(resize,resize);
	M->size = resize;
	return M;
}

//function to free allocated memory from matrix type
void freematrix(int ** m, int size){
	int i;
	for(i=0;i<size;i++)	{
		free(m[i]);
	}
	free(m);
}

//function to free allocated memory from mat type structures
void freestruct2(abmatrix * M){
	freematrix(M->a, M->size);
	freematrix(M->b, M->size);
	free(M);
}

//function to free allocated memory from matS type structures
void freestruct(onematrix * M){
	freematrix(M->m, M->size);
	free(M);
}





int ** addmatrix(int ** a, int ** b, int size){
	int i,j;
	int ** c = creatematrix(size,size);
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			c[i][j]=a[i][j]+b[i][j];
		}
	}
	return c;
}

int ** submatrix(int ** a, int ** b, int size){
	int i,j;
	int ** c = creatematrix(size,size);
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			c[i][j]=a[i][j]-b[i][j];
		}
	}
	return c;
}

//Thread Routine that is called recursively
void * remul(void * inputstruct)
{
	int i,j;
	abmatrix * M = (abmatrix *) inputstruct;
	int size = M->size;
	int resize = size/2;
	
	onematrix * c = initstruct(size); //allocating memory to the resultant matrix
	
	if(size==1){
		c->m[0][0] = M->a[0][0] * M->b[0][0] ;
	}
	else{	
		int ** a11 = creatematrix(resize,resize);
		int ** a12 = creatematrix(resize,resize);
		int ** a21 = creatematrix(resize,resize);
		int ** a22 = creatematrix(resize,resize);
		int ** b11 = creatematrix(resize,resize);
		int ** b12 = creatematrix(resize,resize);
		int ** b21 = creatematrix(resize,resize);
		int ** b22 = creatematrix(resize,resize);
		int ** ctemp;
		//initialising the submatrices
		for(i=0;i<resize;i++){
			for(j=0;j<resize;j++){
				a11[i][j] = M->a[i][j];
				b11[i][j] = M->b[i][j];
			}
		}
		int k=0;
		for(i=resize;i<size;i++){
			for(j=0;j<resize;j++){
				a21[k][j] = M->a[i][j];
				b21[k][j] = M->b[i][j];
			}
			k++;
		}
		int l=0;
		for(i=0;i<resize;i++){
			for(j=resize;j<size;j++){
				a12[i][l] = M->a[i][j];
				b12[i][l] = M->b[i][j];
				l++;
			}
			l=0;
		}
		k=0;
		l=0;
		for(i=resize;i<size;i++){
			for(j=resize;j<size;j++){
				a22[k][l] = M->a[i][j];
				b22[k][l] = M->b[i][j];
				l++;
			}
			l=0;
			k++;
		}
		//Creating the structures which will be used to store the matrix pairs following the algorithm
		abmatrix * m1 = initstruct2(resize);//A11, S1
		m1->a = a11;
		m1->b = submatrix(b12, b22, resize);
		abmatrix * m2 = initstruct2(resize);//S2, B22
		m2->a = addmatrix(a11, a12, resize);
		m2->b = b22;
		abmatrix * m3 = initstruct2(resize);//S3, B11
		m3->a = addmatrix(a21, a22, resize);
		m3->b = b11;
		abmatrix * m4 = initstruct2(resize);//A22, S4	
		m4->a = a22;
		m4->b = submatrix(b21, b11, resize);
		abmatrix * m5 = initstruct2(resize);//S5, S6	
		m5->a = addmatrix(a11, a22, resize);
		m5->b = addmatrix(b11, b22, resize);
		abmatrix * m6 = initstruct2(resize);//S7, S8
		m6->a = submatrix(a12, a22, resize);
		m6->b = addmatrix(b21, b22, resize);
		abmatrix * m7 = initstruct2(resize);//S9, S10	
		m7->a = submatrix(a11, a21, resize);
		m7->b = addmatrix(b11, b12, resize);
		
		//Recursive calls 
		onematrix * p1 = remul((void *) m1);
		onematrix * p2 = remul((void *) m2);
		onematrix * p3 = remul((void *) m3);
		onematrix * p4 = remul((void *) m4);
		onematrix * p5 = remul((void *) m5);
		onematrix * p6 = remul((void *) m6);
		onematrix * p7 = remul((void *) m7);
		
		freestruct2(m1);
		freestruct2(m2);
		freestruct2(m3);
		freestruct2(m4);
		freestruct2(m5);
		freestruct2(m6);
		freestruct2(m7);
		
		//After getting the resultant matrices from the recursive calls, 
		//pair matrices for further operations specified by the algorithm
		int ** t1;
		t1 = addmatrix(p5->m, p4->m, resize);
		int ** t2;
		t2 = addmatrix(t1, p6->m, resize);
		ctemp = submatrix(t2, p2->m, resize);
		for(i=0;i<resize;i++){
			for(j=0;j<resize;j++){
				c->m[i][j] = ctemp[i][j];
			}
		}	
		freematrix(t1, resize);
		freematrix(t2, resize);
		freematrix(ctemp, resize);	
		ctemp = addmatrix(p1->m, p2->m, resize);
		l = 0;
		for(i=0;i<resize;i++){
			for(j=resize;j<size;j++){
				c->m[i][j] = ctemp[i][l];
				l++;
			}
			l = 0;
		}
		freematrix(ctemp, resize);
	
	
		ctemp = addmatrix(p3->m, p4->m, resize);
		k = 0;
		for(i=resize;i<size;i++){
			for(j=0;j<resize;j++){
				c->m[i][j] = ctemp[k][j];		
			}
			k++;
		}
		freematrix(ctemp, resize);	
		
		t1 = addmatrix(p5->m, p1->m, resize);
		t2 = submatrix(t1, p3->m, resize);
		ctemp = submatrix(t2, p7->m, resize);
		k = 0;
		l = 0;
		for(i=resize;i<size;i++){
			for(j=resize;j<size;j++){
				c->m[i][j] = ctemp[k][l];
				l++;
			}
			l = 0;
			k++;
		}
		freematrix(t1, resize);
		freematrix(t2, resize);
		freematrix(ctemp, resize);
		freestruct(p1);
		freestruct(p2);
		freestruct(p3);
		freestruct(p4);
		freestruct(p5);
		freestruct(p6);
		freestruct(p7);		
	}
	return (void *) c;
}
		
//Wrapper function for the routine function,
//Here is where the threads are created and initialised		
onematrix * multiply(abmatrix * ABmatrix)
{
	int i,j;
	int size=ABmatrix->size;
	int resize=size/2;
	pthread_t pt[7];
	
	void * stat1;
	void * stat2;
	void * stat3;
	void * stat4;
	void * stat5;
	void * stat6;
	void * stat7;
	
	onematrix * p1;
	onematrix * p2;
	onematrix * p3;
	onematrix * p4;
	onematrix * p5;
	onematrix * p6;
	onematrix * p7;
	onematrix * c = initstruct(size);
	
	int ** a11 = creatematrix(resize,resize);
	int ** a12 = creatematrix(resize,resize);
	int ** a21 = creatematrix(resize,resize);
	int ** a22 = creatematrix(resize,resize);
	int ** b11 = creatematrix(resize,resize);
	int ** b12 = creatematrix(resize,resize);
	int ** b21 = creatematrix(resize,resize);
	int ** b22 = creatematrix(resize,resize);
	int ** ctemp;
	
		for(i=0;i<resize;i++){
			for(j=0;j<resize;j++){
				a11[i][j] = ABmatrix->a[i][j];
				b11[i][j] = ABmatrix->b[i][j];
			}
		}
		int k=0;
		for(i=resize;i<size;i++){
			for(j=0;j<resize;j++){
				a21[k][j] = ABmatrix->a[i][j];
				b21[k][j] = ABmatrix->b[i][j];
			}
			k++;
		}
		int l=0;
		for(i=0;i<resize;i++){
			for(j=resize;j<size;j++){
				a12[i][l] = ABmatrix->a[i][j];
				b12[i][l] = ABmatrix->b[i][j];
				l++;
			}
			l=0;
		}
		k=0;
		l=0;
		for(i=resize;i<size;i++){
			for(j=resize;j<size;j++){
				a22[k][l] = ABmatrix->a[i][j];
				b22[k][l] = ABmatrix->b[i][j];
				l++;
			}
			l=0;
			k++;
		}
	
	//Splitting the input matrices into their submatrices and pairing them accordingly
	abmatrix * m1 = initstruct2(resize);//A11, S1
	m1->a = a11;
	m1->b = submatrix(b12, b22, resize);
	
	abmatrix * m2 = initstruct2(resize);//S2, B22
	m2->a = addmatrix(a11, a12, resize);
	m2->b = b22;
	
	abmatrix * m3 = initstruct2(resize);//S3, B11
	m3->a = addmatrix(a21, a22, resize);
	m3->b = b11;
	
	abmatrix * m4 = initstruct2(resize);//A22, S4	
	m4->a = a22;
	m4->b = submatrix(b21, b11, resize);
	
	abmatrix * m5 = initstruct2(resize);//S5, S6	
	m5->a = addmatrix(a11, a22, resize);
	m5->b = addmatrix(b11, b22, resize);
	
	abmatrix * m6 = initstruct2(resize);//S7, S8
	m6->a = submatrix(a12, a22, resize);
	m6->b = addmatrix(b21, b22, resize);
	
	abmatrix * m7 = initstruct2(resize);//S9, S10	
	m7->a = submatrix(a11, a21, resize);
	m7->b = addmatrix(b11, b12, resize);
	
	//Creating the threads and passing the matrix pairs to the thread routine
	pthread_create(&pt[0], NULL, remul, (void*) m1);
	pthread_create(&pt[1], NULL, remul, (void*) m2);
	pthread_create(&pt[2], NULL, remul, (void*) m3);
	pthread_create(&pt[3], NULL, remul, (void*) m4);
	pthread_create(&pt[4], NULL, remul, (void*) m5);
	pthread_create(&pt[5], NULL, remul, (void*) m6);
	pthread_create(&pt[6], NULL, remul, (void*) m7);
	
	//Collecting the threads and storing the return values in the second parameter passed
	pthread_join(pt[0], &stat1);
	
	p1 = (onematrix *) stat1; //A11, B11
	pthread_join(pt[1], &stat2);
	
	p2 = (onematrix *) stat2;//A12, B21
	pthread_join(pt[2], &stat3);
	
	p3 = (onematrix *) stat3;//A11, B12
	pthread_join(pt[3], &stat4);
	
	p4 = (onematrix *) stat4;//A12, B22
	pthread_join(pt[4], &stat5);
	
	p5 = (onematrix *) stat5;//A21, B11
	pthread_join(pt[5], &stat6);
	
	p6 = (onematrix *) stat6;//A22, B21
	pthread_join(pt[6], &stat7);
	
	p7 = (onematrix *) stat7;//A21, B12
	
	freestruct2(m1);
	freestruct2(m2);
	freestruct2(m3);
	freestruct2(m4);
	freestruct2(m5);
	freestruct2(m6);
	freestruct2(m7);
	
	//Resultant matrices are paired accordingly for further operations as specified by the algorithm
	int ** t1;
	t1 = addmatrix(p5->m, p4->m, resize);
	int ** t2;
	t2 = addmatrix(t1, p6->m, resize);
	ctemp = submatrix(t2, p2->m, resize);
	for(i=0;i<resize;i++){
			for(j=0;j<resize;j++){
				c->m[i][j] = ctemp[i][j];
		}
	}	
	freematrix(t1, resize);
	freematrix(t2, resize);
	freematrix(ctemp, resize);	
	ctemp = addmatrix(p1->m, p2->m, resize);
	l = 0;
	for(i=0;i<resize;i++){
			for(j=resize;j<size;j++){
				c->m[i][j] = ctemp[i][l];
			l++;
		}
		l = 0;
	}
	freematrix(ctemp, resize);
	
	
	ctemp = addmatrix(p3->m, p4->m, resize);
	k = 0;
	for(i=resize;i<size;i++){
			for(j=0;j<resize;j++){
				c->m[i][j] = ctemp[k][j];		
		}
		k++;
	}
	freematrix(ctemp, resize);	
		
	t1 = addmatrix(p5->m, p1->m, resize);
	t2 = submatrix(t1, p3->m, resize);
	ctemp = submatrix(t2, p7->m, resize);
	k = 0;
	l = 0;
	for(i=resize;i<size;i++){
			for(j=resize;j<size;j++){
			c->m[i][j] = ctemp[k][l];
			l++;
		}
		l = 0;
		k++;
	}
	freematrix(t1, resize);
	freematrix(t2, resize);
	freematrix(ctemp, resize);
	freestruct(p1);
	freestruct(p2);
	freestruct(p3);
	freestruct(p4);
	freestruct(p5);
	freestruct(p6);
	freestruct(p7);	
	return c;
}


int main(){	
	int i,j;
	double start, end;
	double time;
	FILE *infile, *outfile;
	int arow,acol,brow,bcol;
	int size;
	abmatrix * ABmatrix = (abmatrix *) malloc(sizeof(abmatrix));
	
	infile = fopen(INPUTFILE, "r");
	if(infile == NULL) {
    	printf("Error in Opening infile");
    	return EXIT_FAILURE;
	}

	fscanf(infile, "%d %d", &arow, &acol);
	if(arow>=acol){
		size=squarematrix(arow);
	}
	else{
		size=squarematrix(acol);
	}
	ABmatrix->size=size;
	ABmatrix->arow=arow;
	ABmatrix->acol=acol;
	ABmatrix->a=creatematrix(size,size);
	for(i=0;i<arow;i++){
		for(j=0;j<acol;j++){
			fscanf(infile, "%d", &ABmatrix->a[i][j]);	
		}			
	}
	
	fscanf(infile, "%d %d", &brow, &bcol);
	ABmatrix->brow=brow;
	ABmatrix->bcol=bcol;	
	ABmatrix->b=creatematrix(size,size);	
	for(i=0;i<brow;i++){
		for(j=0;j<bcol;j++){
			fscanf(infile, "%d", &ABmatrix->b[i][j]);	
		}			
	}

	onematrix * Cmatrix;
	
	start =  omp_get_wtime();
	
	Cmatrix = multiply(ABmatrix);

	end =  omp_get_wtime();
	
	time=((double)end-start);
	
	outfile = fopen("outputR_pthread.txt", "w+");

	printf("Strassen recursive + pthread: %lf\n",time);
	fprintf(outfile, "Strassen recursive + pthread: %lf s\n", time);
	fprintf(outfile, "%d %d \n", arow, bcol);
	
	for(i=0;i<arow;i++) {
		for(j=0;j<bcol;j++){
				fprintf(outfile, "%8.0d ", Cmatrix->m[i][j]);
		}
		fprintf(outfile, "\n");		
	}
	fclose(outfile);
	freestruct2(ABmatrix);
	freestruct(Cmatrix);
	return 0;
}
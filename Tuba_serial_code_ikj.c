/*


 * cache_blocking_and_loop_unrolling_ikj.c
 *
 *  Created on: Nov 29, 2020
 *      Author: tuba

*/



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#define SIZE 8000

double randBtw(double lowerb, double upperb)
{
	double r = upperb - lowerb;
	double d = RAND_MAX / r;
	return lowerb + (rand() / d);
}


void fill1dMat(double *A, int row,int col)
{
	int i,j;
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			//A[i * col +j] = 1+floor(drand48()*1000);
			A[i * col +j] = randBtw(0,1000);
			//A[i * col +j] = (i+1)*(j+2);
		}
	}
}


void print1dMat(double *A, int row, int col)
{
	int i,j;
	for(i = 0; i < row ; i++)
	{
		for(j = 0; j < col; j++)
		{
			printf("%6.1lf\t", A[i * col + j]);
		}
		printf("\n");
	}
}

void transpose(double *A,int row, int col)
{
	int i,j;
	double tmp;
	for(i = 0; i < row; i++)
	{
		for ( j = i+1; j < col; j++)
		{
			tmp = A[i * col + j];
			A[i * col + j]= A[j * row + i];
			A[j * row + i]=tmp;
		}
	}

}

void MatMatMult(double *A, double *B, double *C, int n)
{

	int i,j,k,m,l,p;
	int nb_of_blocks,block_size;
	block_size=8 ;
	nb_of_blocks= n / block_size;
	//transpose(B,n);
	//double r;
	
	for(i=0 ; i < nb_of_blocks; i++)
	{
		for (l = 0 ; l < nb_of_blocks; l++)
		{
			for ( k=0; k < block_size ; k++)
			{
				for(p=0; p < block_size; p++)
				{
					//r= A[i*block_size*n + l*block_size + k*n +p];
					for ( j=0 ; j<nb_of_blocks ; j++)
						{
						for (m=0; m < block_size ; m+=8)
						{
							//printf("%f\n",B[j*block_size*n + l *block_size + m*n +p]);

							//printf("A[ilkp] * B[ljpm] = %f * %f\n", A[i*block_size*n + l*block_size + k*n +p], B[l*block_size*n + j *block_size + p*n +m]);
							//printf("A[ilkp] * B[ljpm] = %f * %f\n\n", A[i*block_size*n + (l)*block_size + k*n +p] , B[l*block_size*n + (j) *block_size + (p)*n + m+1]);
							//printf("A[ilkp] * B[ljpm] = %f * %f\n\n", A[i*block_size*n + (l)*block_size + k*n +p] , B[l*block_size*n + (j) *block_size + (p)*n + m+19]);

							//C[i*block_size*n + j*block_size + k*n + m] += r * B[j*block_size*n + l *block_size + m*n +p];
							C[i*block_size*n + j*block_size + k*n + m] += A[i*block_size*n + l*block_size + k*n +p] *
									B[l*block_size*n + j *block_size + p*n +m];
							C[i*block_size*n + j*block_size + k*n + m+1] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+1];
							C[i*block_size*n + j*block_size + k*n + m+2] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+2];
							C[i*block_size*n + j*block_size + k*n + m+3] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+3];
							C[i*block_size*n + j*block_size + k*n + m+4] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+4];
							C[i*block_size*n + j*block_size + k*n + m+5] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+5];
							C[i*block_size*n + j*block_size + k*n + m+6] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+6];
							C[i*block_size*n + j*block_size + k*n + m+7] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+7];
									/*
							C[i*block_size*n + j*block_size + k*n + m+8] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+8];
							C[i*block_size*n + j*block_size + k*n + m+9] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+9];
							C[i*block_size*n + j*block_size + k*n + m+10] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+10];
							C[i*block_size*n + j*block_size + k*n + m+11] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+11];
							C[i*block_size*n + j*block_size + k*n + m+12] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+12];
							C[i*block_size*n + j*block_size + k*n + m+13] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+13];
							C[i*block_size*n + j*block_size + k*n + m+14] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+14];
							C[i*block_size*n + j*block_size + k*n + m+15] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+15];
							C[i*block_size*n + j*block_size + k*n + m+16] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+16];
							C[i*block_size*n + j*block_size + k*n + m+17] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+17];
							C[i*block_size*n + j*block_size + k*n + m+18] += A[i*block_size*n + (l)*block_size + k*n +p] *
									B[l*block_size*n + (j) *block_size + (p)*n + m+18];
							C[i*block_size*n + j*block_size + k*n + m+19] += A[i*block_size*n + l*block_size + k*n +p] *
									B[l*block_size*n + j *block_size + p*n +m+19];
							*/
							//printf("C[ijkm] = %f \n",C[i*block_size*n + j*block_size + k*n + m+19]);
							//printf("C[ijkm] = %f \n",C[i*block_size*n + j*block_size + k*n + m+1]);
							//printf("C[ijkm] = %f \n",C[i*block_size*n + j*block_size + k*n + m+2]);
							
						}

					}
					

				}
			}
		}
	}
	
									/*
	for (kk=0; kk<n; kk+=block_size)
	{ 
		for (jj=0; jj<n; jj+=block_size) 
		{
			for (i=0; i<n; i++) 
			{
				for (k=kk; k<kk+block_size; k++) 
				{
					for (j=jj; j<jj+block_size;j++)
					{ 
						C[i*(jj+block_size)+j]+=A[i*(kk+block_size)+k]*B[k*(jj+block_size)+j];
					}
				}
			}
		}
	} */
}

int main(int argc, char *argv[])
{
	double *A, *B, *C;
	struct timeval t;
	double time1, time2, wct;

	A = (double*)malloc(SIZE * SIZE * sizeof(double));
	B = (double*)malloc(SIZE * SIZE * sizeof(double));
	C = (double*)malloc(SIZE * SIZE * sizeof(double));

	fill1dMat(A, SIZE,SIZE);
	fill1dMat(B,SIZE,SIZE);

	/*
	printf("A\n");
	print1dMat(A,SIZE,SIZE);
	printf("B\n");
	print1dMat(B,SIZE,SIZE);
	*/

	//print1dMat(A,SIZE,SIZE);
	gettimeofday(&t, NULL);
	time1 = t.tv_sec + 1.0e-6*t.tv_usec;

	MatMatMult(A,B,C,SIZE);
	//transpose(A,SIZE);

	sleep(5);

	gettimeofday(&t, NULL);
	time2 = t.tv_sec + 1.0e-6*t.tv_usec;

	wct = time2 - time1;

	printf("wct = %lf\n", wct);


	//printf("C\n");
	//print1dMat(C,SIZE,SIZE);



	free(B);
	free(C);
	free(A);

	return 0;

}



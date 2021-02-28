/*

 * cache_blocking_and_loop_unrolling_ijk.c
 *
 *  Created on: Nov 28, 2020
 *      Author: tuba

*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>


#define SIZE 2048


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

void MatMatMult(double *A, double *B, double *C, int n){
	int i,j,k,m,l,p;
	double sum;
	int nb_of_blocks,block_size;
	block_size= 8;
	nb_of_blocks= n / block_size;
	transpose(B,n,n);

	for(i=0 ; i < nb_of_blocks; i++)
	{
		for ( j=0 ; j<nb_of_blocks ; j++)
		{
			for ( k=0; k < block_size ; k++)
			{
				for (m=0; m < block_size ; m++)
				{
					sum =0.0;
					for (l = 0 ; l < nb_of_blocks; l++)
					{
						for(p=0; p < block_size; p+=8)
						{
				sum += A[i*block_size*n + l*block_size + k*n +p] *
					B[j*block_size*n + l *block_size + m*n +p];
				sum += A[i*block_size*n + (l)*block_size + k*n +p+1] *
					B[j*block_size*n + (l) *block_size + m*n +p+1];
				sum += A[i*block_size*n + (l)*block_size + k*n +p+2] *
					B[j*block_size*n + (l) *block_size + m*n +p+2];
				sum += A[i*block_size*n + (l)*block_size + k*n +p+3] *
					B[j*block_size*n + (l) *block_size + m*n +p+3];
				sum += A[i*block_size*n + (l)*block_size + k*n +p+4] *
					B[j*block_size*n + (l) *block_size + m*n +p+4];
				sum += A[i*block_size*n + (l)*block_size + k*n +p+5] *
					B[j*block_size*n + (l) *block_size + m*n +p+5];
				sum += A[i*block_size*n + (l)*block_size + k*n +p+6] *
					B[j*block_size*n + (l) *block_size + m*n +p+6];
				sum += A[i*block_size*n + (l)*block_size + k*n +p+7] *
					B[j*block_size*n + (l) *block_size + m*n +p+7];
						}

					}
				C[i*block_size*n + j*block_size + k*n + m] = sum;

				}
			}
		}
	}
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
	
	gettimeofday(&t, NULL);
	time1 = t.tv_sec + 1.0e-6*t.tv_usec;

	MatMatMult(A,B,C,SIZE);
	//transpose(A,SIZE);

	sleep(5);

	gettimeofday(&t, NULL);
	time2 = t.tv_sec + 1.0e-6*t.tv_usec;

	wct = time2 - time1;

	printf("wct = %lf\n", wct);

	//print1dMat(A,SIZE,SIZE);
	//printf("B\n");
	//print1dMat(B,SIZE,SIZE);
	//printf("C\n");
	//print1dMat(C,SIZE,SIZE);



	free(B);
	free(C);
	free(A);

	return 0;

}









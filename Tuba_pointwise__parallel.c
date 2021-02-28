#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define SIZE 2048
//#define seed 5


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
			//A[i * col +j] = (i+1);
		}
	}
}


void MatMatMult_basic(double *A, double *B, double *C, int row1, int col1, int row2, int col2)
{
	int i,j,k;
	double sum;
	for(i = 0; i < row1; i++)
	{
		for (k = 0; k < col1; k++)
		{
			r = A[i * col1 + k];
			for( j = 0; j < col2; j++)
			{
				C[i * col2 + j] += r * B[k * col2 + j];

			}
		}
	}
}

int main(int argc, char ** argv)
{
	int my_rank, np,i,j;
	int block_size, send_count;
	double *B, *A_block, *B_block, *C_block, *C;
	double start, finish;
	
	

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// number of rows to be processed by each slave
	block_size = SIZE / (np); // blk_size = 20:4=5
	int disp[block_size][np]; 
	int rcv_counts[np];
	//int *disp;

	
	A_block = (double*)malloc(block_size*SIZE* sizeof(double));
	B_block = (double*)malloc(block_size* SIZE*sizeof(double));
	C_block = (double*)malloc(block_size* SIZE*sizeof(double));
	B = (double*)calloc(SIZE * SIZE , sizeof(double));
	//A = (double*)malloc(SIZE * SIZE * sizeof(double));
	if (my_rank == 0)
		C = (double*)calloc(SIZE * SIZE , sizeof(double));

	//disp = (int*)calloc(np*np, sizeof(int));

	//MPI_Barrier(MPI_COMM_WORLD);

	fill1dMat(A_block, block_size, SIZE);
	fill1dMat(B_block, block_size, SIZE);


	/*
	for (i=0; i< np; i++)
		rcv_counts[i]=SIZE*block_size;
		*/

	MPI_Barrier(MPI_COMM_WORLD);


	//printf("\nprocessor %d ,A_block\n", my_rank);
	//print1dMat(A_block, block_size, SIZE);

	
	//MPI_Barrier(MPI_COMM_WORLD);
	
	//printf("\nprocessor %d ,B_block\n", my_rank);
	//print1dMat(B_block, block_size, SIZE);

	//MPI_Barrier(MPI_COMM_WORLD);

	

	for (i=0; i< np; i++)
		rcv_counts[i]=1;


	for (j=0; j<block_size; j++)
	{
		for (i=0; i< np; i++)
			{
				//disp[j*np+i]= i * block_size*SIZE+SIZE*j; 
				disp[j][i]= i * block_size*SIZE+SIZE*j;  //0, 16, ...
				//disp[j][i]= i *SIZE*+SIZE*j;  //0, 16, ...

			}
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	/*
	printf("disp\n");

	for (i = 0; i < block_size; i++)
	{
		for (j = 0; j < np; j++)
		{
			printf(" %d\t", disp[i][j]);
		}
		printf("\n");
	}
	*/
	
		//disp[i]= i*block_size;
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	

	MPI_Datatype new_B_block, B_blk;

	int starts[2] = {0,0};
	int subsizes[2] = {1, SIZE};
	int bigsizes[2] = {block_size, SIZE};
	MPI_Type_create_subarray(2, bigsizes, subsizes, starts, MPI_ORDER_C,  MPI_DOUBLE, &new_B_block);
	MPI_Type_create_resized(new_B_block, 0, sizeof(double), &B_blk);
	MPI_Type_commit(&B_blk);
	
	//MPI_Barrier(MPI_COMM_WORLD);
	
	
	
	for (i=0; i< block_size; i++)
		
		//MPI_Allgatherv(&B_block[block_size*i],1, B_blk, B, rcv_counts, disp[i],  B_blk, MPI_COMM_WORLD );
		MPI_Allgatherv(&B_block[SIZE*i],1, B_blk, B, rcv_counts, disp[i],  B_blk, MPI_COMM_WORLD );
		//MPI_Allgatherv(&B_block[SIZE*i],SIZE, MPI_DOUBLE, B, rcv_counts, disp[i], MPI_DOUBLE, MPI_COMM_WORLD );


	/*
	if (my_rank == 0){
	printf("\n processor %d , B\n", my_rank);

	print1dMat(B, SIZE,SIZE);
	
	}
	*/

	MPI_Barrier(MPI_COMM_WORLD);


	
	
	start = MPI_Wtime();
	//transpose(B,SIZE, SIZE);

	/*
	if (my_rank == 0){
	printf("\n processor %d , B_transpose\n", my_rank);

	print1dMat(B, SIZE,SIZE);
	
	}

	*/

	MatMatMult_basic(A_block, B, C_block, block_size, SIZE,  SIZE, SIZE );
	//if (my_rank == 0)
	//{
		//MatMatMult_b_small(A_block, B, C_block, block_size, SIZE, SIZE);
	//}
	
	//MatMatMult(A_block, B, C_block, block_size, SIZE, SIZE);

	//MPI_Barrier(MPI_COMM_WORLD);
	
	
	//printf("\nprocessor %d ,C_block\n", my_rank);
	//print1dMat(C_block,block_size,SIZE);
	
	//MPI_Barrier(MPI_COMM_WORLD);

	for (i=0; i< block_size; i++)
		
		//MPI_Allgatherv(&C_block[block_size*i],1, B_blk, C, rcv_counts, disp[i],  B_blk, MPI_COMM_WORLD );
		MPI_Gatherv(&C_block[SIZE*i],1, B_blk, C, rcv_counts, disp[i],  B_blk,0, MPI_COMM_WORLD );

	//MPI_Barrier(MPI_COMM_WORLD);

	finish = MPI_Wtime();
	if (my_rank==0)
	{
		printf("Proc %d > Elapsed time = %6.4lf seconds\n", my_rank, finish-start);
		free(C);
	}

	MPI_Barrier(MPI_COMM_WORLD);



	/*
	if (my_rank == 0){
	printf("\n processor %d , C\n", my_rank);

	print1dMat(C, SIZE,SIZE);
	
	}
	*/
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Type_free (&B_blk);


	free(A_block);
	//free(A);
	free(B_block);
	free(C_block);
	free(B);
	

	
	MPI_Finalize();
	return 0;
}


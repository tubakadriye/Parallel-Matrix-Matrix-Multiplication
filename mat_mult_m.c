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


void MatMatMult(double *A, double *B, double *C, int row1, int col1, int col2)
{

	int i,j,k,m,l,p;
	double r;
	int nb_of_blocks1,nb_of_blocks2,blk_size;
	blk_size= 8;
	nb_of_blocks1= row1 / blk_size;
	nb_of_blocks2 = col1 /blk_size;
	//transpose(B,n);
	//double r;
	
	
	for(i=0 ; i < nb_of_blocks1; i++)
	{
		for (l = 0 ; l < nb_of_blocks2; l++)
		{
			for ( k=0; k < blk_size ; k++)
			{
				for(p=0; p < blk_size; p++)
				{
					r=A[i*blk_size*col1 + l*blk_size + k*col1 +p]; 
					//r= A[i*block_size*n + l*block_size + k*n +p];
					for ( j=0 ; j<nb_of_blocks2 ; j++)
					{
						for (m=0; m < blk_size ; m+=8)
						{
							//printf("%f\n",B[j*block_size*n + l *block_size + m*n +p]);

							//printf("A[ilkp] * B[ljpm] = %f * %f\n", A[i*block_size*n + l*block_size + k*n +p], B[l*block_size*n + j *block_size + p*n +m]);
							//printf("A[ilkp] * B[ljpm] = %f * %f\n\n", A[i*block_size*n + (l)*block_size + k*n +p] , B[l*block_size*n + (j) *block_size + (p)*n + m+1]);
							//printf("A[ilkp] * B[ljpm] = %f * %f\n\n", A[i*block_size*n + (l)*block_size + k*n +p] , B[l*block_size*n + (j) *block_size + (p)*n + m+19]);

							//C[i*block_size*n + j*block_size + k*n + m] += r * B[j*block_size*n + l *block_size + m*n +p];
							C[i*blk_size*col2 + j*blk_size + k*col2 + m] += r *
									B[j*blk_size*col2 + l *blk_size + m*col2 +p];
							C[i*blk_size*col2 + j*blk_size + k*col2 + m+1] += r *
									B[j*blk_size*col2 + l *blk_size + (m+1)*col2 +p];
							C[i*blk_size*col2 + j*blk_size + k*col2 + m+2] += r *
									B[j*blk_size*col2 + l *blk_size + (m+2)*col2 +p];
							C[i*blk_size*col2 + j*blk_size + k*col2 + m+3] += r *
									B[j*blk_size*col2 + l *blk_size + (m+3)*col2 +p];
							C[i*blk_size*col2 + j*blk_size + k*col2 + m+4] += r *
									B[j*blk_size*col2 + l *blk_size + (m+4)*col2 +p];
							C[i*blk_size*col2 + j*blk_size + k*col2 + m+5] += r *
									B[j*blk_size*col2 + l *blk_size + (m+5)*col2 +p];
							C[i*blk_size*col2 + j*blk_size + k*col2 + m+6] += r *
									B[j*blk_size*col2 + l *blk_size + (m+6)*col2 +p];
							C[i*blk_size*col2 + j*blk_size + k*col2 + m+7] += r *
									B[j*blk_size*col2 + l *blk_size + (m+7)*col2 +p];
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

	//MatMatMult_basic(A_block, B, C_block, block_size, SIZE,  SIZE, SIZE );
	//if (my_rank == 0)
	//{
		//MatMatMult_b_small(A_block, B, C_block, block_size, SIZE, SIZE);
	//}
	
	MatMatMult(A_block, B, C_block, block_size, SIZE, SIZE);

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


/*
Uses Fox's algorithm to multiply two square matrices
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N 256

struct grid_info {
	int numOfProc;  //total number of processes
	int gOrder;     //Order of grid
	int row_no;     //row number
	int col_no;     //column number
	int myId;       //global process id
	int grid_Id;    //pocess id in grid  
	MPI_Comm comm;  //Communicator for entire grid
	MPI_Comm row_comm;  //Communicator for row
	MPI_Comm col_comm; //Communicator for column
	
};

void createGrid(struct grid_info* grid);
void readMatrix(char* filename, double** local, struct grid_info* grid, int n);
void Fox(struct grid_info* grid, int n, double** localA, double** localB, double** localC);
void partMatrixMultiply(double** localA, double** localB, double** localC, int n) 
void printMatrix(char* filename, double**global, double** local, struct grid_info* grid, int n);

int localNum;

int main(int argc, char* argv[]) {
//	int myId;
//	int localNum;
	double** localA;
	double** localB;
	double** localC;
	double** glocalC;
	double start;
	double end;
	double cost;
	struct grid_info grid;
//	MPI_Status status;
	
	char fileA[] = "./matrixA_256x256.txt";
	char fileB[] = "./matrixB_256x256.txt";
	char fileC[] = "./matrixC_256x256.txt";
	
	//Initialize MPI
	MPI_Init(&argc, &argv);

	//Create 2-dim grid
	createGrid(&grid);
	localNum = N/grid.gOrder;
	localA = (double**)malloc(sizeof(double*) * localNum);
	localB = (double**)malloc(sizeof(double*) * localNum);
	localC = (double**)malloc(sizeof(double*) * localNum);
	int i;
	for (i = 0; i < localNum; ++i) {
		localA[i] = (double*)malloc(sizeof(double) * localNum);
		localB[i] = (double*)malloc(sizeof(double) * localNum;
		localC[i] = (double*)malloc(sizeof(double) * localNum);
	}
	
	globalC = (double**)malloc(sizeof(double*) * N);
	for(i = 0; i < N; ++i) {
		globalC[i] = (double*)malloc(sizeof(double) * N);
	}
	
	//Read matrix localA from file
	Read_matrix(fileA, localA, &grid, localNum);
	//Read matrix localB from file
	Read_matrix(fileB, localB, &grid, localNum);
	
	if (grid.myId == 0) {
		start = MPI_Wtime();
	}
	//Use fox to multiply two matrices
	Fox(&grid, localNum, localA, localB, localC);
	if (grid.myId == 0) {
		end = MPI_Wtime();
		cost = end - cost;
		printf("computation cost for %d x %d is %lf \n", N, N, cost);
	}
	//Save result to file
	printMatrix(fileC, globalC, localC, &grid, localNum);

	free(localA);
	free(localB);
	free(localC);
	free(globalC);
	
	MPI_Finalize();
	return 0;
}

//Create 2-dim grid
void createGrid(struct grid_info* grid) {
//	int oldId;
	int dim[2];
	int wrap[2];
	int coords[2];
	int sub_row[2];
	int sub_col[2];
	
	//Set up global grid info
	MPI_Comm_rank(MPI_COMM_WORLD, &(grid->myId));
	MPI_Comm_size(MPI_COMM_WORLD, &(grid->numOfProc));
	
	//Get order of grid
	grid->gOrder = (int) sqrt(grid->numOfProc);
	dim[0] = grid->gOrder;
	dim[1] = grid->gOrder;
	
	//Cyclic shift
	wrap[0] = 1;
	wrap[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, wrap, 1, &(grid->comm));
	MPI_Comm_rank(grid->comm, &(grid->grid_Id));
	MPI_Cart_coords(grid->comm, grid->grid_Id, 2, coords);
	grid->row_no = coords[0];
	grid->col_no = coords[1];
	sub_row[0] = 0; 
	sub_row[1] = 1;
	sub_col[0] = 1; 
	sub_col[1] = 0;
	
	MPI_Cart_sub(grid->comm, sub_row, &(grid->row_comm));
	MPI_Cart_sub(grid->comm, sub_col, &(grid->col_comm));
}

//process 0 read matrix from file and distribute it to all processors including itself
void readMatrix(char* filename, double** local, struct grid_info* grid, int n) {
	int j, k;
	int index[2];
	int dest;
	MPI_Status status;
	
	if(grid->grid_Id == 0) {
		FILE* fp;
		double** temp;
		int jj = 0;
		int kk = 0;
		double value;
		long int offset;
		temp = (double**)malloc(sizeof(double*) * n);
		for(j = 0; j < n; ++j)
			temp[j] = (double*)malloc(sizeof(double) * n);
		fp = fopen(filename, "r");
		
		for(j = 0; j < grid->gOrder; ++j) {
			index[0] = j;
			for(k = 0; k < grid->gOrder; ++k) {
				index[1] = k;
				MPI_Cart_rank(grid->comm, index, &dest);
				if (dest == 0) {
					for(jj = 0; jj < n; ++jj) {
						fseek(fp, jj * N * sizeof(double), SEEK_SET);
						for (kk = 0; kk < n; ++kk) {
							fscanf(fp, "%lf", &value);
							local[jj][kk] = value;
						}
					}
				}
				else {
					offset = (j * N + k) * n * sizeof(double);
					for(jj = 0; jj < n; ++jj) {
						fseek(fp, offset + jj * N * sizeof(double), SEEK_SET);
						for (kk = 0; kk < n; ++kk) {
							fscanf(fp, "%lf", &value);
							local[jj][kk] = value;
						}
					}
					
					MPI_Send(temp, n * n, MPI_DOUBLE, dest, 0, grid->comm);
				}				
			}			
		}
		
		free(temp);
		
		fclose(fp);
		
	}
	else {
		MPI_Recv(local, localNum * localNum, MPI_DOUBLE, 0, 0, grid->comm, &status);		
	}
	
}

//Fox's algorithm to matrix multiplication
void Fox(struct grid_info* grid, int n, double** localA, double** localB, double** localC) {
	int stage,
	int i, j;
	int localSize = n * n;
	int bcast_root;
	MPI_Status status;
	//Calculate addresses for circular shift of B
	int src = (grid->row_no + 1)%(grid->gOrder);
	int dest = (grid->row_no + grid->gOrder - 1)%(grid->gOrder);
	
	
	tempA = (double**)malloc(sizeof(double*) * localNum);
		for(j = 0; j < localNum; ++j)
			tempA[j] = (double*)malloc(sizeof(double) * localNum);
	
	for(stage = 0; stage < grid->gOrder, ++stages) {
		bcast_root = (grid->row_no + stage) % (grid->gOrder);
		if (bcast_root == grid->col_no) {
			MPI_Bcast(localA, localSize, MPI_DOUBLE, bcast_root, grid->row_comm);
			partMatrixMultiply(localA, localB, localC, localNum);
			
		}
		
		else {
			MPI_Bcast(tempA, localSize, MPI_DOUBLE, bcast_root, grid->row_comm);
			partMatrixMultiply(localA, localB, localC, localNum);
		}
		
		MPI_Sendrecv_replace(localB, localSize, MPI_DOUBLE, dest, 0, src, 0, grid->col_comm, &status);
	}
		
}

//Matrix multiply with two local matrices
void partMatrixMultiply(double** localA, double** localB, double** localC, int n) {
	int i, j, k;
	for(i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			localC[i][j] = 0.0;
			for (k = 0; k < n; ++k) {
				localC[i][j] += localA[i][k]*localB[k][j];
			}
		}		
	}
}

//Print result matrix to file
void printMatrix(char* filename, double**global, double** local, struct grid_info* grid, int n) {
	int index[2];
	int source;
	int i,j;
	int ii,jj;
	int rowIdx, colIdx;
	
	MPI_Status status;
	
	if(grid->myId == 0) {
		double** temp;
		temp = (double**)malloc(sizeof(double*) * n);
		for(j = 0; j < n; ++j)
			temp[j] = (double*)malloc(sizeof(double) * n);
		
		for(i = 0; i < grid->gOrder; ++i) {
			index[0] = i;
			for(j = 0; j < grid->gOrder; ++j) {
				index[1] = j;
				MPI_Cart_rank(grid->comm, index, &source);
				if(source == 0) {
					for(ii = 0; ii < n; ++ii) {
						//rowIdx = ii;
						for(jj = 0; jj < n; ++jj) {
							global[ii][jj] = local[ii][jj];
						}						
					}
					
				}
				else {
					MPI_Recv(temp, n * n; MPI_DOUBLE, source, 1, grid->comm, &status);
					for (ii = 0; ii < n; ++ii) {
						rowIdx = i * n + ii;
						for(jj = 0; jj < n; ++jj) {
							colIdx = j * n + jj;
							global[rowIdx][colIdx] = temp[ii][jj];
						}
					}
				}
			}
			
		}
		free(temp);
	}
	else {
		MPI_Send(local, n * n, MPI_DOUBLE, 0, 1, grid->comm);		
	}	
	
	if (grid->myId == 0) {
		FILE* fp;
		fp = fopen(filename, "w");
	
		if (fp == NULL) {
			printf("cannot open result file!\n");
			exit(-1);
		}
		
		for (i = 0; i < N; ++i) {
			for(j = 0; j < N; ++j) {
				fprintf("%lf\t", global[i][j]);				
			}
			fprintf("\n");
		}
		
		fclose(fp);
	}
	
}




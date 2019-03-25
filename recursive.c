#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<math.h>
#include "mpi.h"
//Method 2: Recursive Doubling
//N: 2^8, 2^9, 2^10, 2^11, 2^12

int main(int argc, char* argv[]){
	int k = 8;
	int numOfProc;
	int myId;
	double start_time, end_time, total_time;
	total_time = 0.0;
   
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &numOfProc);
	MPI_Comm_rank (MPI_COMM_WORLD, &myId);
	
//	char resultFile[] = "./result.txt";
//	FILE *fp = NULL;
	int i,j;
	int idx;
	idx = 0;
	int N = pow(2, k);
	double ii = 0.01;
	double h = 2.0/(N+1);
	int local_columns = N/numOfProc;
	gsl_matrix* local_g_cons = gsl_matrix_alloc(1, local_columns);
	gsl_matrix* local_g_coef = gsl_matrix_alloc(1, local_columns);
	gsl_matrix* local_x_cons = gsl_matrix_alloc(1, local_columns);
	gsl_matrix* local_x_coef = gsl_matrix_alloc(1, local_columns);
	
	gsl_matrix* a = gsl_matrix_alloc(1, N);
	gsl_matrix* b = gsl_matrix_alloc(1, N);
	gsl_matrix* c = gsl_matrix_alloc(1, N);
	gsl_matrix* d = gsl_matrix_alloc(1, N);
	gsl_matrix* w = gsl_matrix_alloc(1, N);
	gsl_matrix* g_cons = gsl_matrix_alloc(1, N);
	gsl_matrix* g_coef = gsl_matrix_alloc(1, N);
	gsl_matrix* x_cons = gsl_matrix_alloc(1, N);
	gsl_matrix* x_coef = gsl_matrix_alloc(1, N);
//Set matrix of a = -k/(h^2)
	for(i = 0; i < N; i++){
		if(i == 0){
			gsl_matrix_set(a, 0, i, 0);
		}
		else{
			gsl_matrix_set(a, 0, i, -ii/(h*h));
		}
	}
	
	//Set matrix of b = 1 + 2*k/(h^2)
	for(i = 0; i < N; i++){
		gsl_matrix_set(b, 0, i, 1+2*ii/(h*h));
	}
	
	//Set matrix of c = -k/(h^2)
	for(i = 0; i < N; i++){
		if(i == N-1){
			gsl_matrix_set(c, 0, i, 0);
		}
		else{
			gsl_matrix_set(c, 0, i, -ii/(h*h));
		}		
	}
	
	//Set matrix of d = u(i, j-1)
	if (idx == 0) {		
		for(i = 0; i < N; i++){
				gsl_matrix_set(d, 0, i, (i+1)*h*(2 - (i+1)*h));
		}
		
	}
	
	//w = ci/(bi-aiwi-1)
	for(i = 0; i < N; i++){
		if(i == 0){
			gsl_matrix_set(w, 0, 0, gsl_matrix_get(c, 0, 0)/gsl_matrix_get(b, 0, 0));
		}
		else{
			gsl_matrix_set(w, 0, i, gsl_matrix_get(c, 0, i)/(gsl_matrix_get(b, 0, i)-gsl_matrix_get(a, 0, i)*gsl_matrix_get(w, 0, i-1)));
		}
	}
	
	//Set matrix of g_cons(0)
	for(i = 0; i < N; i++){
		if(i == 0){
			gsl_matrix_set(g_cons, 0, 0, gsl_matrix_get(d, 0, 0)/gsl_matrix_get(b, 0, 0));
		}
		else{
			gsl_matrix_set(g_cons, 0, i, gsl_matrix_get(d, 0, i)/(gsl_matrix_get(b, 0, i)-gsl_matrix_get(a, 0, i)*gsl_matrix_get(w, 0, i-1)));
		}
	}
	
	//Set matrix of g_coef(0)
	for(i = 0; i < N; i++){
		if(i == 0){
			gsl_matrix_set(g_coef, 0, 0, 0);
		}
		else{
			gsl_matrix_set(g_coef, 0, i, -gsl_matrix_get(a, 0, i)/(gsl_matrix_get(b, 0, i)-gsl_matrix_get(a, 0, i)*gsl_matrix_get(w, 0, i-1)));
		}
	}
	
//	double tmp = 0.0;
	
	if(myId == 0){
		 start_time = MPI_Wtime();
	/*	 fp = fopen(resultFile, "w");
		 if (fp == NULL) {
			 printf("cannot open file\n");
			 exit(-1);			 
		 }
		 for (i = 0; i < N; ++i) {
			 tmp = gsl_matrix_get(d, 0, i);
			 fprintf(fp, "%lf\n", tmp);
		 }
		 fprintf(fp, "\n");,*/
		 
	}
	
	for (idx = 1; idx < 10; ++idx) {
	int start_Index, recursive_interval; 

	for(i = 0; i < k; i++){
		start_Index = pow(2, i);
		recursive_interval = pow(2, i); 
		MPI_Scatter((*g_cons).data, local_columns, MPI_DOUBLE, (*local_g_cons).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter((*g_coef).data, local_columns, MPI_DOUBLE, (*local_g_coef).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		MPI_Bcast((*g_cons).data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*g_coef).data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//In the section where the start_Index located		
		if(myId == start_Index/local_columns){
			for(j = start_Index - myId*local_columns; j < local_columns; j++){	
				//d = g_cons    a = g_coef				
				//dj=ajdj-1+dj
				gsl_matrix_set(local_g_cons, 0, j, gsl_matrix_get(local_g_coef, 0, j)*gsl_matrix_get(g_cons, 0, myId*local_columns+j-recursive_interval)+gsl_matrix_get(local_g_cons, 0, j));
				//aj=aj*aj-1
				gsl_matrix_set(local_g_coef, 0, j, gsl_matrix_get(local_g_coef, 0, j)*gsl_matrix_get(g_coef, 0, myId*local_columns+j-recursive_interval));
			}
		}
		//Beyond the section that the start_Index located, for loop all the elements in the section
		else if(myId > start_Index/local_columns){
			for(j = 0 ; j < local_columns; j++){	
				//d = g_cons    a = g_coef				
				//dj=ajdj-1+dj
				gsl_matrix_set(local_g_cons, 0, j, gsl_matrix_get(local_g_coef, 0, j)*gsl_matrix_get(g_cons, 0, myId*local_columns+j-recursive_interval)+gsl_matrix_get(local_g_cons, 0, j));
				//aj=aj*aj-1
				gsl_matrix_set(local_g_coef, 0, j, gsl_matrix_get(local_g_coef, 0, j)*gsl_matrix_get(g_coef, 0, myId*local_columns+j-recursive_interval));
			}
		}	
		//Before the section that the start_Index located: do nothing
		else{

		} 
		MPI_Gather((*local_g_cons).data, local_columns, MPI_DOUBLE, (*g_cons).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather((*local_g_coef).data, local_columns, MPI_DOUBLE, (*g_coef).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
//**********************************************************************	
	//Now we know g, next we compute Ux = g
	gsl_matrix_memcpy(x_cons, g_cons);
	gsl_matrix_memcpy(x_coef, w);
	gsl_matrix_scale(x_coef, -1);
	if(myId == 0 && k < 5){
	//	print_matrix("x_cons(0): ", x_cons);
		//print_matrix("x_coef(0): ", x_coef);
	}
	for(i = 0; i < k; i++){		
		start_Index = N - 1- pow(2, i);
		recursive_interval = pow(2, i); 
		//In the section that the start_Index located
		MPI_Scatter((*x_cons).data, local_columns, MPI_DOUBLE, (*local_x_cons).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter((*x_coef).data, local_columns, MPI_DOUBLE, (*local_x_coef).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*x_cons).data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*x_coef).data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if(myId == start_Index/local_columns){			
			for( j = 0; j < start_Index - myId*local_columns + 1; j++){
				//dj=ajdj+1+dj
				gsl_matrix_set(local_x_cons, 0, j, gsl_matrix_get(local_x_coef, 0, j)*gsl_matrix_get(x_cons, 0, myId*local_columns+j+recursive_interval)+gsl_matrix_get(local_x_cons, 0, j));
				//aj=aj*aj+1
				gsl_matrix_set(local_x_coef, 0, j, gsl_matrix_get(local_x_coef, 0, j)*gsl_matrix_get(x_coef, 0, myId*local_columns+j+recursive_interval));
			}
		}
		//Before the section that the start_Index located
		else if(myId < start_Index/local_columns){
			for(j = 0; j < local_columns; j++){	
				//d = x_cons    a = x_coef				
				//dj=ajdj+1+dj
				gsl_matrix_set(local_x_cons, 0, j, gsl_matrix_get(local_x_coef, 0, j)*gsl_matrix_get(x_cons, 0, myId*local_columns+j+recursive_interval)+gsl_matrix_get(local_x_cons, 0, j));
				//aj=aj*aj+1
				gsl_matrix_set(local_x_coef, 0, j, gsl_matrix_get(local_x_coef, 0, j)*gsl_matrix_get(x_coef, 0, myId*local_columns+j+recursive_interval));
			}
		}
		else{

		}
		MPI_Gather((*local_x_cons).data, local_columns, MPI_DOUBLE, (*x_cons).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather((*local_x_coef).data, local_columns, MPI_DOUBLE, (*x_coef).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} 
	
 	if(myId == 0 && idx == 9){
		end_time = MPI_Wtime();
    	total_time = (end_time - start_time);
	} 	
	
/*	gsl_matrix_memcpy(d, x_cons);
	for(i = 0; i < N; ++i) {
		for (i = 0; i < N; ++i) {
			 tmp = gsl_matrix_get(d, 0, i);
			 fprintf(fp, "%lf\t", tmp);
		 }
		 fprintf(fp, "\n");
	}  */
	
	}
//	end_time = MPI_Wtime();
 //  	total_time = end_time - start_time;
//	double error = gsl_matrix_get(b, 0, 0)*gsl_matrix_get(x_cons, 0, 0)+ gsl_matrix_get(c, 0, 0)*gsl_matrix_get(x_cons, 0, 1) - gsl_matrix_get(d, 0, 0);
	if(myId == 0){
		printf("Taking tota_time = %f s\n", total_time);
		//fclose(fp);
	}
	MPI_Finalize();
	return 0;
}
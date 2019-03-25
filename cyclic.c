#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<stdlib.h>
#include<math.h>
#include "mpi.h"
//**********************************************
//Method 1: Cyclic Reduction
//**********************************************
//N: 2^8, 2^9, 2^10, 2^11, 2^12

int main(int argc, char* argv[]){
	int numOfProc;
	int myId;
	double start_time, end_time, total_time;
	total_time = 0.0;
    
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &numOfProc);
	MPI_Comm_rank (MPI_COMM_WORLD, &myId);
	
//	char resultFile[] = "./result.txt"; 
	
//	FILE *fp = NULL;
	int idx = 0;
	int i,j;
	int k = 8;
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
//	gsl_matrix* u = gsl_matrix_alloc(N, 10);
	gsl_matrix* p = gsl_matrix_alloc(2*(k+1), N);
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
	
	//Set matrix of g_cons(0) = di/(bi-aiwi-1)
	for(i = 0; i < N; i++){
		if(i == 0){
			gsl_matrix_set(g_cons, 0, 0, gsl_matrix_get(d, 0, 0)/gsl_matrix_get(b, 0, 0));
		}
		else{
			gsl_matrix_set(g_cons, 0, i, gsl_matrix_get(d, 0, i)/(gsl_matrix_get(b, 0, i)-gsl_matrix_get(a, 0, i)*gsl_matrix_get(w, 0, i-1)));
		}
	}
	
	//Set matrix of g_coef(0) = -ai/(bi-aiwi-1)
	for(i = 0; i < N; i++){
		if(i == 0){
			gsl_matrix_set(g_coef, 0, 0, 0);
		}
		else{
			gsl_matrix_set(g_coef, 0, i, -gsl_matrix_get(a, 0, i)/(gsl_matrix_get(b, 0, i)-gsl_matrix_get(a, 0, i)*gsl_matrix_get(w, 0, i-1)));
		}
	}
	
/*	double tmp = 0.0;
	
	if (myId == 0) {
		fp = fopen(resultFile, "w");
		if (fp == NULL) {
			printf("can not open file\n");
			exit(-1);
		}
		for (i = 0; i < N; ++i){
			tmp = gsl_matrix_get(d, 0, i);
			fprintf(fp, "%lf\t", tmp);
		}
		
		fprintf(fp, "\n");
		
	}  */
	
	for (idx = 1; idx < 10; ++idx) {
		
	//Lg=d
	//cyclic reduction for computing g
	if(myId == 0){
		start_time = MPI_Wtime();
	}
	//Set matrix of p of Lu=g
	for(i = 0; i < N; i++){
		gsl_matrix_set(p, 0, i, gsl_matrix_get(g_cons, 0, i));
	}
	for(i = 0; i < N; i++){
		gsl_matrix_set(p, 1, i, gsl_matrix_get(g_coef, 0, i));
	}
	
	int start_Index, local_index, recursive_interval; 
	int e;
	//Lg=d
	//Cyclic reduction for computing g
	for(i = 0; i < k; i++){
		start_Index = pow(2, i)*2-1;
		recursive_interval = pow(2, i+1); 
		MPI_Scatter((*g_cons).data, local_columns, MPI_DOUBLE, (*local_g_cons).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter((*g_coef).data, local_columns, MPI_DOUBLE, (*local_g_coef).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		MPI_Bcast((*g_cons).data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*g_coef).data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//Every processor check if it has the term: (1,3,5,7...)(3,7,11...)
		for(j = 0; start_Index + recursive_interval*j < N; j++){
			if(myId ==  (start_Index+recursive_interval*j)/local_columns){
				local_index = start_Index+recursive_interval*j-myId*local_columns;
				//dj=ajdj-1+dj
				gsl_matrix_set(local_g_cons, 0, local_index, gsl_matrix_get(local_g_coef, 0, local_index)*gsl_matrix_get(g_cons, 0, start_Index+recursive_interval*j-recursive_interval/2)+gsl_matrix_get(local_g_cons, 0, local_index));
				//aj=aj*aj-1
				gsl_matrix_set(local_g_coef, 0, local_index, gsl_matrix_get(local_g_coef, 0, local_index)*gsl_matrix_get(g_coef, 0, start_Index+recursive_interval*j-recursive_interval/2));
			}
			MPI_Gather((*local_g_cons).data, local_columns, MPI_DOUBLE, (*g_cons).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather((*local_g_coef).data, local_columns, MPI_DOUBLE, (*g_coef).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			if(myId == 0){
				for(e = 0; e < N; e++){
					gsl_matrix_set(p, 2*i+2, e, gsl_matrix_get(g_cons, 0, e));
				}
				for(e = 0; e < N; e++){
					gsl_matrix_set(p, 2*i+3, e, gsl_matrix_get(g_coef, 0, e));
				}				
			}
		}
	}

	gsl_matrix* g_value = gsl_matrix_alloc(1, N); 
	gsl_matrix_memcpy(g_value, g_cons);
	
	//Back substitution for incomplete g_value
	if (myId == 0){
		for(i = k-2; i >= 0; i--){
			for(j = pow(2,i) + pow(2,i+1) - 1; j < N; j = j + pow(2, i+1)){
				gsl_matrix_set(g_value, 0, j, (gsl_matrix_get(g_value, 0, j + pow(2,i))-gsl_matrix_get(p, 2*i, j+pow(2, i)))/gsl_matrix_get(p, 2*i+1, j+pow(2, i)));
			}
		}
	}

	//Now we know g, next we compute Ux = g
	//cyclic reduction for computing x (xi=gi-wi*x(i+1))
	gsl_matrix_memcpy(x_cons, g_value);
	gsl_matrix_memcpy(x_coef, w);
	gsl_matrix_scale(x_coef, -1);

	//Set matrix of p of Ux=g
	for(i = 0; i < N; i++){
		gsl_matrix_set(p, 0, i, gsl_matrix_get(x_cons, 0, i));
	}
	for(i = 0; i < N; i++){
		gsl_matrix_set(p, 1, i, gsl_matrix_get(x_coef, 0, i));
	}
	for(i = 0; i < k; i++){		
		start_Index = N - 1- pow(2, i);
		recursive_interval = pow(2, i+1); 
		//In the section that the start_Index located
		MPI_Scatter((*x_cons).data, local_columns, MPI_DOUBLE, (*local_x_cons).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter((*x_coef).data, local_columns, MPI_DOUBLE, (*local_x_coef).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*x_cons).data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*x_coef).data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(j = 0; j < k+1; j++){
			int index = recursive_interval*j;
			int local_index;
			if(myId == index/local_columns && index+recursive_interval/2 < N){
				local_index = index - myId*local_columns;
				//dj=ajdj+1+dj
				gsl_matrix_set(local_x_cons, 0, local_index, gsl_matrix_get(local_x_coef, 0, local_index)*gsl_matrix_get(x_cons, 0, index+recursive_interval/2)+gsl_matrix_get(local_x_cons, 0, local_index));
				//aj=aj*aj+1
				gsl_matrix_set(local_x_coef, 0, local_index, gsl_matrix_get(local_x_coef, 0, local_index)*gsl_matrix_get(x_coef, 0, index+recursive_interval/2));	
			}
		}
		MPI_Gather((*local_x_cons).data, local_columns, MPI_DOUBLE, (*x_cons).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather((*local_x_coef).data, local_columns, MPI_DOUBLE, (*x_coef).data, local_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if(myId == 0){
			for(e = 0; e < N; e++){
				gsl_matrix_set(p, 2*i+2, e, gsl_matrix_get(x_cons, 0, e));
			}
			for(e = 0; e < N; e++){
				gsl_matrix_set(p, 2*i+3, e, gsl_matrix_get(x_coef, 0, e));
			}				
		}
	} 	
	//Back substitution for incomplete g_value
	gsl_matrix* x_value = gsl_matrix_alloc(1, N); 
	gsl_matrix_memcpy(x_value, x_cons);
	if (myId == 0){
		for(i = k-2; i >= 0; i--){
			for(j = N - pow(2,i) - pow(2,i+1); j > 0; j = j - pow(2, i+1)){
				//gsl_matrix_set(x_value, 0, j, (gsl_matrix_get(x_value, 0, j - pow(2,i))-gsl_matrix_get(p, 2*i, j-pow(2,i)))/gsl_matrix_get(p, 2*i+1, j-pow(2,i)));
				gsl_matrix_set(x_value, 0, j, (gsl_matrix_get(d, 0, j - pow(2,i))-gsl_matrix_get(p, 2*i, j-pow(2,i)))/gsl_matrix_get(p, 2*i+1, j-pow(2,i)));
			}
		}
		
		end_time = MPI_Wtime();
    	total_time = (end_time - start_time) + total_time;
		double tm = 0.0;
		for(i = 0; i < N; ++i) {
			if (FP_NAN == isnan(gsl_matrix_get(x_value, 0, i))) {
				tm = i*h*(2-i*h);
				gsl_matrix_set(x_value, 0, i, tm);
			}			
		}
		
		gsl_matrix_memcpy(d, x_value);
		
/*		for (i = 0; i < N; ++i){
			tmp = gsl_matrix_get(d, 0, i);
			fprintf(fp, "%lf\t", tmp);
		}
		
		fprintf(fp, "\n");*/
		
	}  
	
	}
	
	if(myId == 0){
	//	fclose(fp);
		printf("Taking tota_time = %f\n s", total_time);
	}
	MPI_Finalize();
	return 0;
}

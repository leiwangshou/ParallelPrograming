#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define NUM_BLOCKS 1
#define BLOCK_SIZE 256
#define NUM_MEM 32768

//stage one
__global__ void prefixOne(int *in_array, int *out_array, int unsize, int size)
{	
	int tid = threadIdx.x;
	for(int j = 0; j < unsize; j++){
		if (j == 0){
			out_array[tid * unsize] = in_array[tid * unsize];
		} else {
			for(int k = 0; k <= j; k++) {
				out_array[tid * unsize + j] += in_array[tid * unsize + k];
			}
		}
	}
}

//stage two
__global__ void prefixTwo(int *in_array, int unsize, int maxid, int idx)
{
	int tid = threadIdx.x;
	if (tid <= maxid) {
		int maxstep = unsize * (int)(powf(2, idx - 1));
		for(int j = 0; j < maxstep; j++) {
			int startIdx = unsize * (int)(powf(2, idx - 1)) * (1 + 2 * tid);
			in_array[startIdx + j] = in_array[startIdx - 1] + in_array[startIdx + j];
		}		
	}
	
	
}

void prefixsum(int blocks, int threads, int steps, int *array_h, int size)
{
	int *array_d;
	int *tmp_one;
	int unsize = size/(blocks * threads);
	
	
	dim3 dim_grid(blocks, 1, 1);
	dim3 dim_block(threads, 1, 1);

	// allocate tmp_d
	cudaMalloc((void **)&tmp_one, size * sizeof(int));
	//cudaMalloc((void **)&out_array_d, blocks * sizeof(int));
	cudaMalloc((void **)&array_d, size * sizeof(int));
	//copy data from host to device
	cudaMemcpy(array_d, array_h, size * sizeof(int),
		   cudaMemcpyHostToDevice);
	
	cudaMemset(tmp_one, 0, size * sizeof(int));
	//do stage 1
	prefixOne<<<dim_grid, dim_block>>> (array_d, tmp_one, unsize, size);
	
	if (steps !=0) {
		int maxtid = 0;
		//do stage 2
		for (int i = 1; i <= steps; i++) {
			maxtid = (int)pow(2, steps-i) - 1;
			prefixTwo<<<dim_grid, dim_block>>>(tmp_one, unsize, maxtid, i);
		}
	}
	
		cudaMemcpy(array_h, tmp_one, size * sizeof(int), cudaMemcpyDeviceToHost);
	
	cudaFree(array_d);
	cudaFree(tmp_one);
}



void prepare_numbers(int **array, int count)
{
	int *numbers = (int *)malloc(count * sizeof(int));

	// load array
	for (int i = 0; i < count; i++) {
		numbers[i] = 1;
	}

	*array = numbers;
}

void print_array(int *array, int count)
{
	for (int i = 0; i < count; i++) {
		printf("%d\t", array[i]);
	}
	printf("\n");
}

int main()
{
	int blocks, threads, max, stepTwo;
	int *array;
    float calTime;
	cudaEvent_t start, end;
	cudaEventCreate(&start);
	cudaEventCreate(&end);

	blocks = NUM_BLOCKS;
	threads = BLOCK_SIZE;
	stepTwo = 8;
	max = NUM_MEM;

	// pre-init numbers
	array = NULL;
	prepare_numbers(&array, max);

	cudaEventRecord(start, 0);
	prefixsum(blocks, threads, stepTwo, array, max);
	cudaEventRecord(end, 0);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&calTime, start, end);

	// print array
//	print_array(array, max);
	printf("the elapsed time with %d threads is %.10f\n", threads, calTime);

	free(array);

	return 0;
}

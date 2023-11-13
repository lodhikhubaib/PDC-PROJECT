#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include "omp.h"

void quickSort_parallel(long long int* array, int lenArray, int numThreads);
void quickSort_parallel_internal(long long int* array, int left, int right, int cutoff);
void quickSort(long long int* arr, int left, int right);


void quickSort_parallel(long long int* array, int lenArray, int numThreads){

	int cutoff = 1000;

	#pragma omp parallel num_threads(numThreads)
	{	
		#pragma omp single nowait
		{
			quickSort_parallel_internal(array, 0, lenArray-1, cutoff);	
		}
	}	

}



void quickSort_parallel_internal(long long int* array, int left, int right, int cutoff) 
{
	
	int i = left, j = right;
	int tmp;
	int pivot = array[(left + right) / 2];

	
	{
	  	/* PARTITION PART */
		while (i <= j) {
			while (array[i] < pivot)
				i++;
			while (array[j] > pivot)
				j--;
			if (i <= j) {
				tmp = array[i];
				array[i] = array[j];
				array[j] = tmp;
				i++;
				j--;
			}
		}

	}


	if ( ((right-left)<cutoff) ){
		if (left < j){ quickSort(array, left, j); }			
		if (i < right){ quickSort(array, i, right); }

	}else{
		#pragma omp task 	
		{ quickSort_parallel_internal(array, left, j, cutoff); }
		#pragma omp task 	
		{ quickSort_parallel_internal(array, i, right, cutoff); }		
	}

}


int lenArr = 100000000;
int numthreads = 4;



void quickSort(long long int* arr, int left, int right) 
{
	int i = left, j = right;
	int tmp;
	int pivot = arr[(left + right) / 2];

  	/* PARTITION PART */
	while (i <= j) {
		while (arr[i] < pivot)
			i++;
		while (arr[j] > pivot)
			j--;
		if (i <= j) {
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	}

	/* RECURSION PART */
	if (left < j){ quickSort(arr, left, j);  }
	if (i< right){ quickSort(arr, i, right); }
}





// Comparator used in qsort()
int cmpfunc (const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}
long long int* arr1;

void check(long long int *arr,int l)
{
for(int i=0;i<l-1;i++)
{
if(arr[i]>arr[i+1])
{printf("Array sorted incorrectly\n");free(arr);exit(0);}
}

}

int main(int argc, char** argv){
	int minMum = 1;
	int maxNum = lenArr;
	numthreads=atoi(argv[1]);
	if(argc==3)
	lenArr=atoi(argv[2]);
	double startTime, stopTime;
	for(int l=1024;l<=lenArr;l*=2){
	arr1 = (long long int*) malloc(l*sizeof(long long int));
	
	// Initialise the array with random numbers
	int i;
	srand(5); // seed
	//printf("Initializing the arrays with random numbers...\n");
	for (i=0; i<l; i++){
		// RAND_MAX = 2147483647 = INT_MAX 
		// printf("RAND_MAX %u \n", RAND_MAX );
		arr1[i] = minMum+(rand()%maxNum);
		//printf("%d \n", arr1[i] ); 
	}

	startTime = omp_get_wtime();
	quickSort_parallel(arr1, l, numthreads);
	stopTime = omp_get_wtime();
	check(arr1,l);
	printf("%10d %10f\n", l,(stopTime-startTime));

	free(arr1);}
	return 0;
}

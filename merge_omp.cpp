#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#define THREADS 4
#define MAX 100000000

unsigned int rand_interval(unsigned int min, unsigned int max)
{
    int r;
    const unsigned int range = 1 + max - min;
    const unsigned int buckets = RAND_MAX / range;
    const unsigned int limit = buckets * range;

    do
    {
        r = rand();
    } 
    while (r >= limit);

    return min + (r / buckets);
}

void fillupRandomly (int *m, int size, unsigned int min, unsigned int max){
    for (int i = 0; i < size; i++)
    m[i] = rand_interval(min, max);
} 

void mergeSortAux(int *X, int n, int *tmp) {
   int i = 0;
   int j = n/2;
   int ti = 0;

   while (i<n/2 && j<n) {
      if (X[i] < X[j]) {
         tmp[ti] = X[i];
         ti++; i++;
      } else {
         tmp[ti] = X[j];
         ti++; j++;
      }
   }
   while (i<n/2) { /* finish up lower half */
      tmp[ti] = X[i];
      ti++; i++;
   }
   while (j<n) { /* finish up upper half */
      tmp[ti] = X[j];
      ti++; j++;
   }
   memcpy(X, tmp, n*sizeof(int));
} 

void mergeSort(int *X, int n, int *tmp,int l)
{
   if (n < 2) return;

   #pragma omp task shared(X) if (n > l)
   mergeSort(X, n/2, tmp,l);

   #pragma omp task shared(X) if (n > l)
   mergeSort(X+(n/2), n-(n/2), tmp + n/2,l);

   #pragma omp taskwait
   mergeSortAux(X, n, tmp);
}

void init(int *a, int size){
   for(int i = 0; i < size; i++)
       a[i] = 0;
}

void printArray(int *a, int size){
   for(int i = 0; i < size; i++)
       printf("%d ", a[i]);
   printf("\n");
}

int isSorted(int *a, int size){
   for(int i = 0; i < size - 1; i++)
      if(a[i] > a[i + 1])
	{printf("Array not sorted correctly\n");
	exit(0);
	return 0;}
   return 1;
}

int main(int argc, char *argv[]) {
        srand(123456);
        int N  = (argc > 1) ? atoi(argv[1]) : MAX;
        int print = (argc > 2) ? atoi(argv[2]) : 0;
        int numThreads = (argc > 3) ? atoi(argv[3]) : THREADS;
        omp_set_dynamic(0);              /** Explicitly disable dynamic teams **/
        omp_set_num_threads(numThreads); /** Use N threads for all parallel regions **/
	int *X,*tmp;	
	printf("%10s | %10s\n\n","SIZE","TIME");
	for(int l=1024;l<MAX;l*=2)
	{
	X = malloc(l * sizeof(int));
        tmp = malloc(l * sizeof(int));

        fillupRandomly (X, l, 1, 10000);

        clock_t t1 = clock();
        #pragma omp parallel
        {
            
            mergeSort(X, l, tmp,l);
        }   
        clock_t t2 = clock();
	isSorted(X, l);
	printf("%10d %10f\n",l,(t2-t1)/(double)CLOCKS_PER_SEC);
        //printf("Time: %f (s) \n",(end-begin)/(double)CLOCKS_PER_SEC);
        free(X);
        free(tmp);
	}
        return (EXIT_SUCCESS);
}

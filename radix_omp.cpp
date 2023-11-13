#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#define BASE_BITS 8
#define L 0.6
#define BASE (1 << BASE_BITS)
#define MASK (BASE-1)
#define DIGITS(v, shift) (((v) >> shift) & MASK)
#define LENGTH 100000000


void omp_lsd_radix_sort(size_t n, unsigned data[n]) {
    unsigned * buffer = malloc(n*sizeof(unsigned));
    int total_digits = sizeof(unsigned)*8;
 
    //Each thread use local_bucket to move data
    size_t i;
    for(int shift = 0; shift < total_digits; shift+=BASE_BITS) {
        size_t bucket[BASE] = {0};
 
        size_t local_bucket[BASE] = {0}; // size needed in each bucket/thread
        //1st pass, scan whole and check the count
        #pragma omp parallel firstprivate(local_bucket)
        {
            #pragma omp for schedule(static) nowait
            for(i = 0; i < n; i++){
                local_bucket[DIGITS(data[i], shift)]++;
            }
            #pragma omp critical
            for(i = 0; i < BASE; i++) {
                bucket[i] += local_bucket[i];
            }
            #pragma omp barrier
            #pragma omp single
            for (i = 1; i < BASE; i++) {
                bucket[i] += bucket[i - 1];
            }
            int nthreads = omp_get_num_threads();
            int tid = omp_get_thread_num();
            for(int cur_t = nthreads - 1; cur_t >= 0; cur_t--) {
                if(cur_t == tid) {
                    for(i = 0; i < BASE; i++) {
                        bucket[i] -= local_bucket[i];
                        local_bucket[i] = bucket[i];
                    }
                } else { //just do barrier
                    #pragma omp barrier
                }
 
            }
            #pragma omp for schedule(static)
            for(i = 0; i < n; i++) { //note here the end condition
                buffer[local_bucket[DIGITS(data[i], shift)]++] = data[i];
            }
        }
        //now move data
        unsigned* tmp = data;
        data = buffer;
        buffer = tmp;
    }
    free(buffer);
}

void isSorted(unsigned *arr,int n)
{
	for(int i=0;i<n-1;i++)
	{
		if(arr[i]>arr[i+1])
		{
			printf("Incorrectly sorted array\n");free(arr);exit(0);
		}
	}
}

int main()
{
	
	clock_t t1,t2;
	printf("%10s | %10s\n\n","SIZE","TIME");
	for(int l=1024;l<=LENGTH;l*=2){
	unsigned * arr = malloc(l*sizeof(unsigned));
	t1=clock();
	omp_lsd_radix_sort(l, arr);
	t2=clock();
	isSorted(arr,l);
	free(arr);
	printf("%10d %10f\n",l,L*(double)(t2-t1)/CLOCKS_PER_SEC);}

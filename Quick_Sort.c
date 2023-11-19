#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#define THREADS 4
#define MAX 1000000

unsigned int rand_interval(unsigned long long int *seed, unsigned int min, unsigned int max) {
    int r;
    const unsigned int range = 1 + max - min;
    const unsigned int buckets = RAND_MAX / range;
    const unsigned int limit = buckets * range;

    do {
        r = rand_r((unsigned int *)seed);  // Casting seed to unsigned int *
    } while (r >= limit);

    return min + (r / buckets);
}

void fillupRandomly(unsigned long long int *m, int size, unsigned int min, unsigned int max, unsigned long long int *seed) {
    for (int i = 0; i < size; i++)
        m[i] = rand_interval(seed, min, max);
}

int isSorted(unsigned long long int *a, int size) {
    for (int i = 0; i < size - 1; i++)
        if (a[i] > a[i + 1]) {
            printf("Array not sorted correctly\n");
            exit(0);
            return 0;
        }
    return 1;
}

void printArray(unsigned long long int *arr, int size) {
    printf("Sorted Array: ");
    for (int i = 0; i < size; i++) {
        printf("%lld ", arr[i]);
    }
    printf("\n");
}

void printOutput(int size, double serialTime, double parallelTime, double speedup) {
    printf("%10d | Serial Time: %.4f seconds | Parallel Time: %.4f seconds | Speedup: %.3f\n",
           size, serialTime, parallelTime, speedup);
}

void swap(long long int *a, long long int *b) {
    long long int t = *a;
    *a = *b;
    *b = t;
}

int partition(long long int arr[], int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (arr[j] <= pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quickSort(long long int arr[], int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void parallelquickSort(long long int arr[], int left, int right) {
    if (left < right) {
        int pivot = arr[left];
        int i = left, j = right;

        while (i < j) {
            while (arr[j] > pivot)
                j--;
            while (i < j && arr[i] <= pivot)
                i++;
            if (i < j) {
                swap(&arr[i], &arr[j]);
            }
        }

        arr[left] = arr[i];
        arr[i] = pivot;

        #pragma omp task
        {
            parallelquickSort(arr, left, i - 1);
        }

        #pragma omp task
        {
            parallelquickSort(arr, i + 1, right);
        }
        #pragma omp taskwait
    }
}

int main(int argc, char *argv[]) {
    int N = (argc > 1) ? atoi(argv[1]) : MAX;
    int print = (argc > 2) ? atoi(argv[2]) : 0;
    int numThreads = (argc > 3) ? atoi(argv[3]) : THREADS;

    if (argc < 4) {
        printf("Usage: %s <array_size> <print_option> <num_threads>\n", argv[0]);
        return EXIT_FAILURE;
    }

    omp_set_num_threads(numThreads);
    omp_set_nested(1);
    omp_set_dynamic(0); // Disable dynamic adjustment of threads
    srand(123456);

    printf("-----------------------------------------------------Quick Sort----------------------------------------------\n\n");
    printf("%10s |      %10s            |       %10s           | %10s \n\n", "SIZE", "SERIAL_TIME", "PARALLEL_TIME", "SPEEDUP");

    long long int *X, *arr;

    for (unsigned long int l = 1024; l < N; l *= 2) {

        X = malloc(l * sizeof(long long int));

        unsigned long long int seed = 123456;

        fillupRandomly(X, l, 1, 10000, &seed);

        // serial
        arr = malloc(l * sizeof(long long int));
        for (int i = 0; i < l; i++) {
            arr[i] = X[i];
        }

        double t1_serial = omp_get_wtime();
        for (int i = 0; i < l; i++) {
            quickSort(arr, 0, l - 1);
        }
        double t2_serial = omp_get_wtime();
        isSorted(arr, l);
        double time_serial = t2_serial - t1_serial;

        // PARALLEL
        arr = malloc(l * sizeof(long long int));
        for (int i = 0; i < l; i++) {
            arr[i] = X[i];
        }

        double t1_parallel = omp_get_wtime();
        #pragma omp parallel
        {
            #pragma omp single
            {
                parallelquickSort(arr, 0, l - 1);
            }
        }
        double t2_parallel = omp_get_wtime();
        isSorted(arr, l);
        double time_parallel = t2_parallel - t1_parallel;
        double speedup = time_serial / time_parallel;

        printOutput(l, time_serial, time_parallel, speedup);

        if (print == 1) {
            printArray(arr, l);
        }

        free(X);
        free(arr);
    }

    return EXIT_SUCCESS;
}

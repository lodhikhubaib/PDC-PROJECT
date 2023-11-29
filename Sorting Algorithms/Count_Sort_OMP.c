#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#define THREADS 4
#define MAX 1000000

unsigned int rand_interval(unsigned int *seed, unsigned int min, unsigned int max) {
    int r;
    const unsigned int range = 1 + max - min;
    const unsigned int buckets = RAND_MAX / range;
    const unsigned int limit = buckets * range;

    do {
        r = rand_r(seed);  // Casting seed to unsigned int *
    } while (r >= limit);

    return min + (r / buckets);
}

void fillupRandomly(unsigned int *m, int size, unsigned int min, unsigned int max, unsigned int *seed) {
    for (int i = 0; i < size; i++)
        m[i] = rand_interval(seed, min, max);
}

int isSorted(unsigned int *a, int size) {
    for (int i = 0; size > 1 && i < size - 1; i++)
        if (a[i] > a[i + 1]) {
            printf("Array not sorted correctly\n");
            exit(0);
            return 0;
        }
    return 1;
}

void printArray(unsigned int *arr, int size) {
    printf("Sorted Array: ");
    for (int i = 0; i < size; i++) {
        printf("%u ", arr[i]);
    }
    printf("\n");
}

void printOutput(int size, double serialTime, double parallelTime, double speedup) {
    printf("%10d | Serial Time: %.4f seconds | Parallel Time: %.4f seconds | Speedup: %.3f\n",
           size, serialTime, parallelTime, speedup);
}

void printArrayToFile(unsigned int *arr, int size, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < size; i++) {
        if (fprintf(file, "%d ", arr[i]) < 0) {
            perror("Error writing to file");
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}


void countingSort(unsigned int A[], int n, int num) {
    int *B = (int *)malloc(n * sizeof(unsigned int));
    for(int i = 0; i < n; i++){
        B[i] = 0;
    }
    int i, j, x, count;
    #pragma omp parallel private(i, j, count) num_threads(num)
    {
        #pragma omp for
        for (i = 0; i < n; i++)
        {
            count = 0;
            for (j = 0; j < n; j++)
            {
                if (A[i] > A[j])
                    count++;
            }
            while (B[count] != 0)
                count++;
            B[count] = A[i];
        }
    }
    for(int i = 0; i < n; i++){
        A[i] = B[i];
    }
    free(B);
}

int main(int argc, char *argv[]) {
    int N = (argc > 1) ? atoi(argv[1]) : MAX;
    int print = (argc > 2) ? atoi(argv[2]) : 0;
    int numThreads = (argc > 3) ? atoi(argv[3]) : THREADS;

    if (argc < 4) {
        printf("Usage: %s <array_size> <print_option> <num_threads>\n", argv[0]);
        return EXIT_FAILURE;
    }

    omp_set_nested(1);
    omp_set_dynamic(0); // Disable dynamic adjustment of threads
    srand(123456);

    printf("-----------------------------------------------------Counting Sort----------------------------------------------\n\n");
    printf("%10s |      %10s            |       %10s           | %10s \n\n", "SIZE", "SERIAL_TIME", "PARALLEL_TIME", "SPEEDUP");

    unsigned int *X, *arr;

    for (unsigned long int l = 1024; l < N; l *= 2) {

        X = malloc(l * sizeof(unsigned int));

        unsigned int seed = 123456;

        fillupRandomly(X, l, 1, 10000, &seed);

        // serial
        arr = malloc(l * sizeof(unsigned int));
        for (int i = 0; i < l; i++) {
            arr[i] = X[i];
        }

        double t1_serial = omp_get_wtime();
        countingSort(arr, l, 1);
        double t2_serial = omp_get_wtime();
        isSorted(arr, l);
        //printf("hi");
        double time_serial = t2_serial - t1_serial;

        // PARALLEL
        arr = malloc(l * sizeof(unsigned int));
        for (int i = 0; i < l; i++) {
            arr[i] = X[i];
        }

        double t1_parallel = omp_get_wtime();
        countingSort(arr, l, numThreads);
        double t2_parallel = omp_get_wtime();
        isSorted(arr, l);
        //printf("Hello");
        double time_parallel = t2_parallel - t1_parallel;
        double speedup = time_serial / time_parallel;

        printOutput(l, time_serial, time_parallel, speedup);

        if (print == 1) {
            printArrayToFile(arr, l, "Count_sorted_output.txt");
        }

        free(X);
        free(arr);
    }

    return EXIT_SUCCESS;
}

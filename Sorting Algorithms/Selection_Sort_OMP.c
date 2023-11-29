#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#define THREADS 4
#define MAX 1000000

struct Compare {
    int val;
    int index;
};

// Custom reduction for finding the index of the max element.
#pragma omp declare reduction(maximum : struct Compare : omp_out = omp_in.val > omp_out.val ? omp_in : omp_out)

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

// Function to check if an array is sorted
int isSorted(unsigned int *a, int size) {
    for (int i = 0; size > 1 && i < size - 1; i++)
        if (a[i] < a[i + 1]) {
            printf("Array not sorted correctly %u %u\n", a[i], a[i + 1]);
            return 0;
        }
    return 1;
}
void parallelSelectionSort(int *A, int n, int numThreads) {
    // Iterate through the array
    for (int startpos = 0; startpos < n; startpos++) {
        struct Compare max;
        max.val = A[startpos];
        max.index = startpos;

        // Divide the array into chunks for parallel processing
#pragma omp parallel num_threads(numThreads)
        {
            struct Compare privateMax = max;

            // Determine the range of elements to process in parallel
#pragma omp for
            for (int i = startpos + 1; i < n; ++i) {
                if (A[i] > privateMax.val) {
                    privateMax.val = A[i];
                    privateMax.index = i;
                }
            }

            // Update the shared max with the maximum from this chunk
#pragma omp critical
            {
                if (privateMax.val > max.val) {
                    max = privateMax;
                }
            }
        }

        // Swap the found maximum element with the first element
        int temp = A[startpos];
        A[startpos] = A[max.index];
        A[max.index] = temp;
    }
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

int main(int argc, char *argv[]) {
    int N = (argc > 1) ? atoi(argv[1]) : MAX;
    int print = (argc > 2) ? atoi(argv[2]) : 0;
    int numThreads = (argc > 3) ? atoi(argv[3]) : THREADS;

    if (argc < 4) {
        printf("Usage: %s <array_size> <print_option> <num_threads>\n", argv[0]);
        return EXIT_FAILURE;
    }

    omp_set_nested(1);
    omp_set_dynamic(1); // Disable dynamic adjustment of threads
    srand(123456);

    printf("-----------------------------------------------------Selection Sort----------------------------------------------\n\n");
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
        parallelSelectionSort(arr, l,1);
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
        parallelSelectionSort(arr, l,numThreads);
        double t2_parallel = omp_get_wtime();
        isSorted(arr, l);
        //printf("Hello");
        double time_parallel = t2_parallel - t1_parallel;
        double speedup = time_serial / time_parallel;

        printOutput(l, time_serial, time_parallel, speedup);

        if (print == 1) {
            printArrayToFile(X, l, "Selection_sorted_output.txt");
        }

        free(X);
        free(arr);
    }

    return EXIT_SUCCESS;
}

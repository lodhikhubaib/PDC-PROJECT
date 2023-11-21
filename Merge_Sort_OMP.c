#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

#define THREADS 4
#define MAX 100000000

unsigned int rand_interval(unsigned int *seed, unsigned int min, unsigned int max) {
    int r;
    const unsigned int range = 1 + max - min;
    const unsigned int buckets = RAND_MAX / range;
    const unsigned int limit = buckets * range;

    do {
        r = rand_r(seed);
    } while (r >= limit);

    return min + (r / buckets);
}

void fillupRandomly(int *m, int size, unsigned int min, unsigned int max, unsigned int *seed) {
    for (int i = 0; i < size; i++)
        m[i] = rand_interval(seed, min, max);
}

void merge(int *X, int left, int mid, int right, int *tmp) {
    int i, j, k;
    int n1 = mid - left + 1;
    int n2 = right - mid;

    int *L = malloc(n1 * sizeof(int));
    int *R = malloc(n2 * sizeof(int));

    for (i = 0; i < n1; i++)
        L[i] = X[left + i];
    for (j = 0; j < n2; j++)
        R[j] = X[mid + 1 + j];

    i = 0;
    j = 0;
    k = left;

    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            X[k] = L[i];
            i++;
        } else {
            X[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        X[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        X[k] = R[j];
        j++;
        k++;
    }

    free(L);
    free(R);
}

void mergeSort(int *X, int n, int *tmp) {
    for (int current_size = 1; current_size < n; current_size = 2 * current_size) {
#pragma omp parallel for schedule(dynamic)//
        for (int left_start = 0; left_start < n; left_start += 2 * current_size) {
            int mid = left_start + current_size - 1;
            int right_end = (left_start + 2 * current_size - 1 < n - 1) ? (left_start + 2 * current_size - 1) : (n - 1);
            merge(X, left_start, mid, right_end, tmp);
        }
    }
}

int isSorted(int *a, int size) {
    for (int i = 0; i < size - 1; i++)
        if (a[i] > a[i + 1]) {
            printf("Array not sorted correctly\n");
            exit(0);
            return 0;
        }
    return 1;
}

void printOutput(int size, double serialTime, double parallelTime, double speedup) {
    printf("%10d | Serial Time: %.4f seconds | Parallel Time: %.4f seconds | Speedup: %.3f\n",
           size, serialTime, parallelTime, speedup);
}

void printArray(int *arr, int size) {
    printf("Sorted Array: ");
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
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
    srand(123456);
    int N = (argc > 1) ? atoi(argv[1]) : MAX;
    int print = (argc > 2) ? atoi(argv[2]) : 0;
    int numThreads = (argc > 3) ? atoi(argv[3]) : THREADS;
    if (argc < 4) {
        printf("Usage: %s <array_size> <print_option> <num_threads>\n", argv[0]);
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);
    print = atoi(argv[2]);
    numThreads = atoi(argv[3]);
    omp_set_dynamic(0);
    omp_set_num_threads(numThreads);

    int *X, *tmp;
    printf("-----------------------------------------------------Merge Sort----------------------------------------------\n\n");
    printf("%10s |      %10s            |       %10s           | %10s \n\n", "SIZE", "SERIAL_TIME", "PARALLEL_TIME" ,"SPEEDUP");
    for (int l = 1024; l < N; l *= 2) {
        X = malloc(l * sizeof(int));
        tmp = malloc(l * sizeof(int));

        unsigned int seed = 123456;

        fillupRandomly(X, l, 1, 10000, &seed);

        // serial
        clock_t t1_serial = clock();
        mergeSort(X, l, tmp);
        clock_t t2_serial = clock();
        isSorted(X, l);
        double time_serial = (t2_serial - t1_serial) / (double)CLOCKS_PER_SEC;

        // PARALLEL
        X = malloc(l * sizeof(int));
        tmp = malloc(l * sizeof(int));

        seed = 123456;

        fillupRandomly(X, l, 1, 10000, &seed);

        clock_t t1_parallel = clock();
#pragma omp parallel
        {
#pragma omp single
            mergeSort(X, l, tmp);
        }
        clock_t t2_parallel = clock();
        isSorted(X, l);
        double time_parallel = (t2_parallel - t1_parallel) / (double)CLOCKS_PER_SEC;
        double speedup = time_serial / time_parallel;

        printOutput(l, time_serial, time_parallel, speedup);

        if(print==1){
            printArrayToFile(X, l, ",Merge_sorted_output.txt");
        }

        free(X);
        free(tmp);
    }

    return EXIT_SUCCESS;
    
}

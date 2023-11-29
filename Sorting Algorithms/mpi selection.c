/*
	Selection Sort for p > 2
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	Syleopoulos Anastasios
*/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h> // sleep

#define _MAXNUMBER 10000 /* maximum array element */

/************** FUNCTIONS ****************/
// fills the array with rand()
void FilltheArray(int *arr, int N) {
    int i;

    for (i = 0; i < N; i++)
        arr[i] = rand() % _MAXNUMBER;
}

// prints the array
void Myprint(int *arr, int N, char *str) {
    int i;

    printf("\n%s\n", str);
    for (i = 0; i < N; i++)
        printf("%d |", arr[i]);
    printf("\n");
}

/***************** MAIN **********************/
int main(int argc, char *argv[]) {
    int myid, numprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int i, j;  // indices for loops
    int temp;   // variable for sorting
    int N, part; // N: size of the array | part: size of each part
    int *index;  // array of indices
    int *A, *As; // A: unsorted array | As: sorted array
    int *Aw;     // array for each part used in selection sort
    double start_time, end_time; // start_time: start time | end_time: end time

    /********** MASTER ***********/
    if (myid == 0) {
        printf("Selection Sort: Parallel program for p > 2\n");
        printf("Enter the size of the array: ");
        scanf("%d", &N);

        A = malloc(N * sizeof(int));
        As = malloc(N * sizeof(int));
        if (A == NULL || As == NULL) {
            printf("\nMemory allocation error!\nEnd of program!\n");
            exit(1);
        }
        FilltheArray(A, N); // fill the array with elements

        part = N / numprocs;  // calculate the part size
        if (N % numprocs != 0) {
            printf("\nThe number of processes (numprocs) is not a multiple of the size (N) of the array!\nEnd of program!\n");
            exit(1);
        }
        printf("The number of processes (numprocs): %d, The size (N) of the array: %d, The part size (part): %d\n", numprocs, N, part);

        start_time = MPI_Wtime(); // start timing

        index = malloc(numprocs * sizeof(int));
        if (index == NULL) {
            printf("\nMemory allocation error!\nEnd of program!\n");
            exit(1);
        }
        for (i = 0; i < numprocs; i++)  // Assign values to the array (index) with the initial position of each node element
            index[i] = i * part;
    }

    MPI_Bcast(&part, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the part size to all processes
    Aw = malloc(part * sizeof(int));                 // Allocate memory for the array used in each process for selection sort
    MPI_Scatter(&A[0], part, MPI_INT, &Aw[0], part, MPI_INT, 0, MPI_COMM_WORLD); // Divide the unsorted array (A) into the arrays Aw

    /************ Selection Sort *************/
    // Sorting of the array in each part
    for (i = 0; i < part; i++) {
        for (j = i + 1; j < part; j++) {
            if (Aw[i] > Aw[j]) {
                temp = Aw[i];
                Aw[i] = Aw[j];
                Aw[j] = temp;
            }
        }
    }

    /****************************************/
    MPI_Gather(&Aw[0], part, MPI_INT, &A[0], part, MPI_INT, 0, MPI_COMM_WORLD); // Gather elements from Aw and store them in array A

    // Merge parts and perform the final sorting of the array
    /********** MASTER ***********/
    if (myid == 0) {
        int tp;  // tp: index of the node
        int check;  // assigned value that is expected to be the smallest
        int done; // check if the smallest has been found (done=1)

        /* 
           For each index of the sorted array (As), find
           the smallest element of the array Aw of each node. 
           Assign the value to the array (As) and increase the index 
           from the node that found the element by one.  
        */ 

        for (i = 0; i < N; i++) {
            tp = 0;
            done = 0;
            while (done != 1) {
                if (index[tp] < (tp + 1) * part) {
                    check = index[tp];
                    for (j = tp + 1; j < numprocs; j++) {
                        if (A[check] > A[index[j]] && index[j] < (j + 1) * part) {
                            check = index[j];
                            tp = j;
                        }
                    }
                    done = 1;
                } else {
                    tp++;
                }
            }
            As[i] = A[check];
            index[tp]++;
        }

        end_time = MPI_Wtime(); // end timing
        printf("\nExecution time: %.16f\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
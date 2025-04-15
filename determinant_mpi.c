#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

// Macros for matrix access using 1D representation
#define MAT(mat, i, j, size) (mat[(i) * (size) + (j)])

// Allocate a 1D matrix (size x size)
long double* allocate_matrix(int size) {
    return (long double*)malloc(size * size * sizeof(long double));
}

// Free allocated matrix memory
void free_matrix(long double* matrix) {
    free(matrix);
}

// Distributed LU decomposition and determinant computation
long double calculate_determinant_distributed(long double* matrix, int size, int rank, int size_comm) {
    long double* L = allocate_matrix(size);
    long double* U = allocate_matrix(size);

    for (int i = 0; i < size * size; i++) {
        L[i] = 0.0L;
        U[i] = 0.0L;
    }

    for (int k = 0; k < size; k++) {
        // Compute U[k][j] row by owning process
        if (k % size_comm == rank) {
            for (int j = k; j < size; j++) {
                long double sum = 0.0L;
                for (int s = 0; s < k; s++)
                    sum += MAT(L, k, s, size) * MAT(U, s, j, size);
                MAT(U, k, j, size) = MAT(matrix, k, j, size) - sum;
            }
        }

        // Broadcast full U[k] row to all processes
        MPI_Bcast(&MAT(U, k, 0, size), size, MPI_LONG_DOUBLE, k % size_comm, MPI_COMM_WORLD);

        // Compute L[i][k] column elements
        for (int i = k; i < size; i++) {
            if (i == k) {
                if (rank == (i % size_comm))
                    MAT(L, i, i, size) = 1.0L;
            } else if (i % size_comm == rank) {
                long double sum = 0.0L;
                for (int s = 0; s < k; s++)
                    sum += MAT(L, i, s, size) * MAT(U, s, k, size);
                if (MAT(U, k, k, size) == 0.0L) return 0.0L; // Singular
                MAT(L, i, k, size) = (MAT(matrix, i, k, size) - sum) / MAT(U, k, k, size);
            }
        }

        // Broadcast L[k+1..n][k] column to all processes in one call
        for (int i = k + 1; i < size; i++) {
            MPI_Bcast(&MAT(L, i, k, size), 1, MPI_LONG_DOUBLE, i % size_comm, MPI_COMM_WORLD);
        }
    }

    // Each process calculates part of the product of U's diagonal
    long double local_det = 1.0L;
    for (int i = 0; i < size; i++) {
        if (i % size_comm == rank)
            local_det *= MAT(U, i, i, size);
    }

    long double global_det;
    MPI_Reduce(&local_det, &global_det, 1, MPI_LONG_DOUBLE, MPI_PROD, 0, MPI_COMM_WORLD);

    free_matrix(L);
    free_matrix(U);

    return global_det;
}

int main(int argc, char* argv[]) {
    int rank, size_comm;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size_comm);

    int matrix_size;
    long double* matrix = NULL;

    if (rank == 0) {
        FILE* file = NULL;
        if (argc > 1) {
            file = fopen(argv[1], "r");
        } else {
            file = fopen("matrix.txt", "r");
        }

        if (!file) {
            fprintf(stderr, "Error: Cannot open matrix file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fscanf(file, "%d", &matrix_size);
        matrix = allocate_matrix(matrix_size);

        srand((unsigned int)time(NULL));  // Seed the random number generator
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                // Generate a random number in the range [-0.1, 0.1]
                MAT(matrix, i, j, matrix_size) = ((rand() % 2001) - 1000) / 20000.0L;
            }
        }


        fclose(file);
    }

    // Broadcast matrix size
    MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate matrix on other ranks
    if (rank != 0)
        matrix = allocate_matrix(matrix_size);

    // Broadcast full matrix row by row
    for (int i = 0; i < matrix_size; i++) {
        MPI_Bcast(&MAT(matrix, i, 0, matrix_size), matrix_size, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Measure time
    double start_time = MPI_Wtime();
    long double determinant = calculate_determinant_distributed(matrix, matrix_size, rank, size_comm);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("\nDeterminant = %.10Le\n", determinant);
        printf("Elapsed Time = %.6f seconds\n", end_time - start_time);
    }

    free_matrix(matrix);
    MPI_Finalize();
    return 0;
}

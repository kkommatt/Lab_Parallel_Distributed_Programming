#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

// Function declarations
long double calculate_determinant(long double** matrix, int size);
long double** allocate_matrix(int size);
void free_matrix(long double** matrix, int size);
int lu_decomposition(long double** matrix, long double** L, long double** U, int size);

int main(int argc, char* argv[]) {
    int matrix_size, num_threads = 16;
    long double** matrix;

    // Set number of OpenMP threads
    if (argc > 2) {
        num_threads = atoi(argv[2]);  // Get number of threads from command line
    }

    omp_set_num_threads(num_threads);  // Apply OpenMP thread count

    // Read matrix from file or generate it randomly
    if (argc > 1) {
        // Read from file
        FILE* file = fopen(argv[1], "r");
        if (!file) {
            fprintf(stderr, "Error: Unable to open file %s\n", argv[1]);
            return 1;
        }

        // Read matrix size and then elements
        fscanf(file, "%d", &matrix_size);
        matrix = allocate_matrix(matrix_size);
        for (int i = 0; i < matrix_size; i++)
            for (int j = 0; j < matrix_size; j++)
                fscanf(file, "%Lf", &matrix[i][j]);

        fclose(file);
    } else {
        // Random matrix generation
        printf("Generating a random matrix. Enter size: ");
        scanf("%d", &matrix_size);

        matrix = allocate_matrix(matrix_size);
        srand((unsigned int)time(NULL));  // Seed random generator

        // Fill matrix with values in range [-1.0, 1.0]
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                int val = rand() % 2001 - 1000;     // Random int between -1000 and 1000
                matrix[i][j] = val / 1000.0L;        // Convert to long double in [-1.0, 1.0]
            }
        }
    }

    // Measure execution time of determinant calculation
    double start = omp_get_wtime();
    long double det = calculate_determinant(matrix, matrix_size);
    double end = omp_get_wtime();

    // Output results
    printf("\nDet = %.10Le\n", det);
    printf("Delta-time: %.6f seconds\n", end - start);

    // Free matrix memory
    free_matrix(matrix, matrix_size);
    return 0;
}

// Allocates a 2D square matrix of size `size × size`
long double** allocate_matrix(int size) {
    long double** mat = (long double**)malloc(size * sizeof(long double*));
    for (int i = 0; i < size; i++)
        mat[i] = (long double*)malloc(size * sizeof(long double));
    return mat;
}

// Frees a previously allocated 2D matrix
void free_matrix(long double** matrix, int size) {
    for (int i = 0; i < size; i++)
        free(matrix[i]);
    free(matrix);
}

// Performs LU Decomposition: matrix = L * U
// Returns 1 if successful, 0 if matrix is singular
int lu_decomposition(long double** matrix, long double** L, long double** U, int size) {
    for (int i = 0; i < size; i++) {
        // Construct upper triangular matrix U
        #pragma omp parallel for  // Parallelize row-wise operation
        for (int k = i; k < size; k++) {
            long double sum = 0.0L;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = matrix[i][k] - sum;
        }

        // Construct lower triangular matrix L
        for (int k = i; k < size; k++) {
            if (i == k)
                L[i][i] = 1.0L;  // Diagonal elements of L are 1
            else {
                long double sum = 0.0L;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);
                if (U[i][i] == 0.0L) return 0;  // Division by zero → singular matrix
                L[k][i] = (matrix[k][i] - sum) / U[i][i];
            }
        }
    }
    return 1;
}

// Computes determinant using LU decomposition
long double calculate_determinant(long double** matrix, int size) {
    long double** L = allocate_matrix(size);
    long double** U = allocate_matrix(size);
    int success = lu_decomposition(matrix, L, U, size);
    long double det = 1.0L;

    if (!success) {
        fprintf(stderr, "Matrix is singular, determinant is 0.\n");
        det = 0.0L;
    } else {
        // Determinant is product of U[i][i]
        for (int i = 0; i < size; i++)
            det *= U[i][i];
    }

    free_matrix(L, size);
    free_matrix(U, size);
    return det;
}

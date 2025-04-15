#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Function prototypes
long double** allocate_matrix(int size);                      // Allocate memory for a square matrix
void free_matrix(long double** matrix, int size);             // Free allocated memory for a matrix
int lu_decomposition(long double** A, long double** L, long double** U, int size);  // Perform LU decomposition
long double calculate_determinant(long double** matrix, int size);                  // Calculate matrix determinant

int main(int argc, char* argv[]) {
    int size;
    long double** matrix;

    // If a filename is passed as argument, read the matrix from the file
    if (argc > 1) {
        FILE* file = fopen(argv[1], "r");
        if (!file) {
            fprintf(stderr, "Error: Unable to open file %s\n", argv[1]);
            return 1;
        }

        fscanf(file, "%d", &size);                      // Read matrix size
        matrix = allocate_matrix(size);                 // Allocate memory for matrix

        // Read matrix values from file
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                fscanf(file, "%Lf", &matrix[i][j]);

        fclose(file);                                   // Close the file
    } else {
        // If no file is provided, ask user for matrix size and generate random values
        printf("Generating a random matrix. Enter size: ");
        scanf("%d", &size);
        matrix = allocate_matrix(size);                 // Allocate memory for matrix

        srand((unsigned int)time(NULL));                // Seed the random number generator

        // Fill matrix with random values scaled to [-0.01, 0.01]
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int val = rand() % 2001 - 1000;          // Random integer in range [-1000, 1000]
                matrix[i][j] = val / 100000.0L;          // Convert to long double in [-0.01, 0.01]
            }
        }
    }

    // Time measurement begins
    clock_t start = clock();

    // Calculate determinant using LU decomposition
    long double det = calculate_determinant(matrix, size);

    // Time measurement ends
    clock_t end = clock();

    // Print results
    printf("Determinant of the matrix: %.30Le\n", det);
    printf("Elapsed time: %.6lf seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    free_matrix(matrix, size);                          // Free allocated memory
    return 0;
}

// Allocate memory for a size x size matrix of long doubles
long double** allocate_matrix(int size) {
    long double** mat = malloc(size * sizeof(long double*));
    for (int i = 0; i < size; i++)
        mat[i] = malloc(size * sizeof(long double));
    return mat;
}

// Free memory allocated for a size x size matrix
void free_matrix(long double** matrix, int size) {
    for (int i = 0; i < size; i++)
        free(matrix[i]);
    free(matrix);
}

// Perform LU decomposition of matrix A into L (lower) and U (upper)
// Returns 1 if successful, 0 if matrix is singular or nearly singular
int lu_decomposition(long double** A, long double** L, long double** U, int size) {
    for (int i = 0; i < size; i++) {
        // Compute upper triangular matrix U
        for (int k = i; k < size; k++) {
            long double sum = 0.0L;
            for (int j = 0; j < i; j++)
                sum += L[i][j] * U[j][k];               // Sum of products of L row and U column
            U[i][k] = A[i][k] - sum;                    // Subtract sum from original matrix entry
        }

        // Compute lower triangular matrix L
        for (int k = i; k < size; k++) {
            if (i == k)
                L[i][i] = 1.0L;                          // Diagonal of L is 1
            else {
                long double sum = 0.0L;
                for (int j = 0; j < i; j++)
                    sum += L[k][j] * U[j][i];           // Sum for L[k][i] calculation
                if (fabsl(U[i][i]) < 1e-20L) return 0;  // Check for near-zero pivot
                L[k][i] = (A[k][i] - sum) / U[i][i];     // Compute L[k][i]
            }
        }
    }
    return 1;
}

// Calculate the determinant of a matrix using LU decomposition
long double calculate_determinant(long double** matrix, int size) {
    long double** L = allocate_matrix(size);
    long double** U = allocate_matrix(size);

    int success = lu_decomposition(matrix, L, U, size); // Decompose matrix
    long double det = 1.0L;

    if (!success) {
        printf("Matrix is singular or near-singular.\n");
        det = 0.0L;
    } else {
        // Determinant is the product of U's diagonal elements
        for (int i = 0; i < size; i++)
            det *= U[i][i];
    }

    free_matrix(L, size);       // Free L matrix
    free_matrix(U, size);       // Free U matrix
    return det;
}

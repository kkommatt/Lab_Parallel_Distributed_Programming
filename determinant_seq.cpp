#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

// Function to calculate the determinant of a matrix
double CalculateDeterminant(vector<vector<double>>& matrix, int matrix_size);

// Function to get the minor of the matrix (excluding the row and column specified)
vector<vector<double>> GetMinor(const vector<vector<double>>& matrix, int row, int col, int matrix_size);

int main(int argc, char* argv[]) {
    int matrix_size;
    vector<vector<double>> matrix;
    if (argc > 1) {
        ifstream file(argv[1]);
        if (file.is_open()) {
            file >> matrix_size;
            matrix.resize(matrix_size, vector<double>(matrix_size));
            for (int i = 0; i < matrix_size; i++) {
                for (int j = 0; j < matrix_size; j++) {
                    file >> matrix[i][j];
                }
            }
            file.close();
        } else {
            cerr << "Error: Unable to open file " << argv[1] << endl;
            return 1;
        }
    } else {
        // Generate a random matrix
        cout << "Generating a random matrix. Enter size: ";
        cin >> matrix_size;
        matrix.resize(matrix_size, vector<double>(matrix_size));
        srand(time(0));
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                matrix[i][j] = rand() % 10 + 1;
            }
        }
    }

    cout << "Input Matrix:" << endl;
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    clock_t start_time = clock();
    double determinant = CalculateDeterminant(matrix, matrix_size);
    clock_t end_time = clock();

    double elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Determinant of the matrix: " << determinant << endl;
    cout << "Elapsed time: " << elapsed_time << " seconds" << endl;

    return 0;
}

// Recursive function to compute determinant
double CalculateDeterminant(vector<vector<double>>& matrix, int matrix_size) {
    if (matrix_size == 1) return matrix[0][0];
    if (matrix_size == 2) return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];

    double determinant = 0;
    for (int col = 0; col < matrix_size; col++) {
        vector<vector<double>> minor = GetMinor(matrix, 0, col, matrix_size);
        determinant += pow(-1.0, col) * matrix[0][col] * CalculateDeterminant(minor, matrix_size - 1);
        cout << "det_i = " << determinant << endl;
    }
    return determinant;
}

// Function to extract minor matrix
vector<vector<double>> GetMinor(const vector<vector<double>>& matrix, int row, int col, int matrix_size) {
    vector<vector<double>> minor(matrix_size - 1, vector<double>(matrix_size - 1));
    int minor_row = 0;
    for (int i = 0; i < matrix_size; i++) {
        if (i == row) continue;
        int minor_col = 0;
        for (int j = 0; j < matrix_size; j++) {
            if (j == col) continue;
            minor[minor_row][minor_col++] = matrix[i][j];
        }
        minor_row++;
    }
    return minor;
}

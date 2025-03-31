#include <iostream>
#include <cmath>
#include <mpi.h>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>

#define MASTER_PROCESS 0
#define MSG_FROM_MASTER 1
#define MSG_FROM_WORKER 2

using namespace std;

double CalculateDeterminant(vector<vector<double> > &matrix, int start_row, int end_row, int matrix_size) {
    if (matrix_size == 1) return matrix[0][0];
    if (matrix_size == 2) return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];

    double determinant = 0;
    for (int column = start_row; column < end_row; column++) {
        vector<vector<double> > sub_matrix(matrix_size - 1, vector<double>(matrix_size - 1));
        for (int row = 1; row < matrix_size; row++) {
            int sub_matrix_column = 0;
            for (int col = 0; col < matrix_size; col++) {
                if (col == column) continue;
                sub_matrix[row - 1][sub_matrix_column++] = matrix[row][col];
            }
        }
        determinant += pow(-1.0, column) * matrix[0][column] * CalculateDeterminant(
            sub_matrix, 0, matrix_size - 1, matrix_size - 1);
    }
    return determinant;
}

int main(int argc, char *argv[]) {
    int total_processes, process_rank, total_workers, source, destination, message_type, extra_rows, row_offset;
    double segment_determinant, start_time, end_time, matrix_read_start_time, matrix_read_end_time;
    int matrix_size;

    char hostname[MPI_MAX_PROCESSOR_NAME];
    int hostname_len;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    MPI_Get_processor_name(hostname, &hostname_len);

    start_time = MPI_Wtime();
    total_workers = total_processes - 1;

    if (process_rank == MASTER_PROCESS) {
        cout << "Determinant calculation" << endl;
        cout << endl << endl;

        matrix_read_start_time = MPI_Wtime();
        vector<double> matrix_data;

        if (argc > 1) {
            ifstream matrix_file(argv[1]);
            matrix_file >> matrix_size;
            matrix_data.resize(matrix_size * matrix_size);
            for (int i = 0; i < matrix_size * matrix_size; i++) {
                matrix_file >> matrix_data[i];
            }
            matrix_file.close();
        } else {
            cout << "Enter matrix size: ";
            cin >> matrix_size;
            matrix_data.resize(matrix_size * matrix_size);
            srand(time(0));
            for (int i = 0; i < matrix_size * matrix_size; i++) {
                matrix_data[i] = rand() % 10 - 5;
            }
        }
        matrix_read_end_time = MPI_Wtime();

        cout << "Input matrix:" << endl;
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                cout << matrix_data[i * matrix_size + j] << " ";
            }
            cout << endl;
        }
        cout << endl << endl;

        float temp = static_cast<float>(matrix_size) / total_processes + 0.5;
        row_offset = static_cast<int>(temp);
        extra_rows = matrix_size % total_processes;

        message_type = MSG_FROM_MASTER;
        for (destination = 1; destination <= total_workers; destination++) {
            MPI_Send(&matrix_size, 1, MPI_INT, destination, message_type, MPI_COMM_WORLD);
            MPI_Send(matrix_data.data(), matrix_size * matrix_size, MPI_DOUBLE, destination, message_type,
                     MPI_COMM_WORLD);
        }

        vector<vector<double> > matrix(matrix_size, vector<double>(matrix_size));
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                matrix[i][j] = matrix_data[i * matrix_size + j];
            }
        }

        double total_determinant = CalculateDeterminant(matrix, 0, row_offset, matrix_size);
        cout << "det_i = " << total_determinant << endl;

        message_type = MSG_FROM_WORKER;
        for (int i = 1; i <= total_workers; i++) {
            MPI_Recv(&segment_determinant, 1, MPI_DOUBLE, i, message_type, MPI_COMM_WORLD, &status);
            total_determinant += segment_determinant;
        }

        end_time = MPI_Wtime();
        cout << endl << endl;
        cout << "Delta-time: " << (end_time - start_time) - (matrix_read_end_time - matrix_read_start_time) <<
                " seconds" << endl;
        cout << endl << endl;

        cout << "Det = " << total_determinant << endl;
    }

    if (process_rank > 0) {
        message_type = MSG_FROM_MASTER;
        MPI_Recv(&matrix_size, 1, MPI_INT, MASTER_PROCESS, message_type, MPI_COMM_WORLD, &status);
        vector<double> matrix_data(matrix_size * matrix_size);
        MPI_Recv(matrix_data.data(), matrix_size * matrix_size, MPI_DOUBLE, MASTER_PROCESS, message_type,
                 MPI_COMM_WORLD, &status);

        float temp = static_cast<float>(matrix_size) / total_processes + 0.5;
        row_offset = static_cast<int>(temp);

        int start_row = process_rank * row_offset;
        int end_row = (process_rank == total_workers) ? matrix_size : start_row + row_offset;

        vector<vector<double> > matrix(matrix_size, vector<double>(matrix_size));
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                matrix[i][j] = matrix_data[i * matrix_size + j];
            }
        }

        segment_determinant = CalculateDeterminant(matrix, start_row, end_row, matrix_size);
        cout << "det_i = " << segment_determinant << endl;

        message_type = MSG_FROM_WORKER;
        MPI_Send(&segment_determinant, 1, MPI_DOUBLE, MASTER_PROCESS, message_type, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

# Determinant Calculation Project

This repository contains implementations of determinant calculation for square matrices using different approaches:
- **Sequential Computation**
- **Parallel Processing with OpenMP**
- **Distributed Computation with MPI**

## Features
- **Reads matrix from file or generates a random matrix**
- **Supports different computation models for performance comparison**
- **Measures execution time for performance analysis**

## Requirements
### Sequential Version
- C++ compiler (GCC, Clang, MSVC)

### OpenMP Version
- Compiler supporting OpenMP (e.g., `g++ -fopenmp`)

### MPI Version
- MPI library (e.g., OpenMPI or MPICH)
- `mpic++` for compilation

## Installation
```sh
git clone https://github.com/your-repo/determinant-calculation.git  
cd determinant-calculation  
```

## Usage
### Sequential
```sh
g++ -o sequential determinant_sequential.cpp  
./sequential matrix.txt  
```

### OpenMP
```sh
g++ -fopenmp -o openmp determinant_openmp.cpp  
./openmp matrix.txt [num_threads]  
```

### MPI
```sh
mpic++ -o mpi determinant_mpi.cpp  
mpirun -np [num_processes] ./mpi matrix.txt  
```

## Performance Considerations
- OpenMP improves performance by leveraging multi-threading on shared memory architectures.
- MPI is ideal for distributed systems where processes communicate over a network.
- Recursive determinant calculation has exponential complexity; alternative algorithms like LU decomposition are recommended for large matrices.

## License
This project is licensed under the MIT License.
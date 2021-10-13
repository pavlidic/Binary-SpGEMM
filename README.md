# Binary-SpGEMM
An MPI+OpenMP hybrid implementation of Gustavson's Sparse Matrix Multiplication for Boolean CSR matrices

Only `final/` and `Matlab/` are intended for use, `old/` has some old code that is there for archiving purposes.

### Compilation ###
Go to `final/` and do `make`

### Usage ###
For performance testing the algorithm:

`mpirun/srun -n number_of_tasks SpGEMM_mpi_omp path_to_matrix block_size number_of_threads times_to_run`

`number_of_tasks*path_to_matrix*block_size` needs to equal the input's size `n`.

It prints in order:

`tasks, threads, total_cpus/threads, blocksize, matrix path, input size n, input nnz, output nnz, mean time, median time, fastest time`

For testing the correctness of the algorithm, use `Matlab/test_SpGEMM(n,d)` where `n` is the size of the `n*n` random input sparse matrix and `d` is aproximately
the non zero elements per row that the matrix will have.

`old/` has the code of the serial implementation + previous implementations that used blocking but were much slower.

SHELL := /bin/bash

CC=mpicc

CFLAGS = -O3 -fopenmp -march=native -w
EXTRALIBS = mmio.o utils.o coo2csc.o
BINS = SpGEMM_mpi_omp

all: $(BINS)

test: SpGEMM_mpi_omp_validity
	mpirun -n 4 SpGEMM_mpi_omp_validity ../Matlab/validity_test.mtx 6250 2

%.o: %.c
	$(CC) $(CFLAGS) -c $^

SpGEMM_mpi_omp: SpGEMM_mpi_omp.o $(EXTRALIBS)
	$(CC) $(CFLAGS) -o $@ $^
	strip $@

SpGEMM_mpi_omp_validity: SpGEMM_mpi_omp_validity.o $(EXTRALIBS)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f *.o *.gcda

clean_all:
	rm -f *.o $(BINS) *.gcda SpGEMM_mpi_omp_validity
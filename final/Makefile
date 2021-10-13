SHELL := /bin/bash

CC=mpicc

CFLAGS = -O3 -fopenmp -march=native -w
EXTRALIBS = mmio.o utils.o coo2csc.o
BINS = SpGEMM_mpi_omp

all: $(BINS)

%.o: %.c
	$(CC) $(CFLAGS) -c $^

SpGEMM_mpi_omp: SpGEMM_mpi_omp.o $(EXTRALIBS)
	$(CC) $(CFLAGS) -o $@ $^
	strip $@

clean:
	rm -f *.o *.gcda

clean_all:
	rm -f *.o $(BINS) *.gcda
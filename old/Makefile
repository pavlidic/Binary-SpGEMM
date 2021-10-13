SHELL := /bin/bash

CC=gcc
MPICC=mpicc

CFLAGS = -O3
OMPFLAGS = -fopenmp
EXTRALIBS = mmio.o utils.o ins.o coo2csc.o
BINS = BSpGEMM SpGEMM_omp

all: $(BINS)

mpi: SpGEMM_mpi_omp.c mmio.c utils.c ins.c coo2csc.c
	mpicc -fopenmp -O3 -o mpitest SpGEMM_mpi_omp.c mmio.c utils.c ins.c coo2csc.c -w -march=native

%.o: %.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c $^

BSpGEMM: BSpGEMM.o $(EXTRALIBS)
	$(CC) $(CFLAGS) -o $@ $^

SpGEMM_omp: SpGEMM_omp.o $(EXTRALIBS)
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $^

clean:
	rm -f *.o *.gcda

clean_all:
	rm -f *.o $(BINS) *.gcda
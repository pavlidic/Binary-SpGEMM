#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>
#include <stdint.h>

#define tic tictoc(0)
#define toc tictoc(1)

void printCSR(int* row, int* col, int n, int m, int b);


void readCOO(const char *mat, uint32_t **row, uint32_t **col,uint32_t *M, uint32_t *N, uint32_t *nnz);

struct timespec diff(struct timespec start, struct timespec end);

double timeConv(struct timespec result);

double tictoc(int mode);

// function to swap elements
void swap(int *a, int *b);

// function to find the partition position
int partition(int array[], int low, int high);

void quickSort(int array[], int low, int high);

// function to swap elements
void swapD(double *a, double *b);

// function to find the partition position
int partitionD(double array[], int low, int high);

void quickSortD(double array[], int low, int high);



#endif


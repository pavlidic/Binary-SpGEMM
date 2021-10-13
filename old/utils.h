#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>
#include <stdint.h>

#define tic tictoc(0)
#define toc tictoc(1)

typedef struct CSRbCSR{
    int  n,m;
    int* LLcol;
    int* LLrow;
    int* Bcol;
    int* Brow;
    int  nnz;
    int  nnzb;
    int  b;
    int  bpr;
    int  bpc;
}CSRbCSR;

int rndCSR(int** row, int** col, int n, int m, int nnz_per_row);

void full2CSR(const bool *C, const int n, const int m, int *Ccol, int *Crow);

bool full2CSR_re(const bool *C, const int n, const int m, int **Ccol, int *Crow, int *Csize);

void full2CSR_T(const bool *C, const int n, const int m, int *Ccol, int *Crow);

void printCSR(int* row, int* col, int n, int m, int b);

void csr2bcsr(int*  csr_row, int*  csr_col, int n, int m, int nnz, int b, CSRbCSR* A, int print);

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

//merges 2 arrays into a third, no doubles, returns n3
int merge(int* arr1, int* arr2, int* arr3, int n1, int n2);

// ORs 2 sparse matrices(CSR)
void SpM_OR(int* Acol, int* Arow, int n, int* Bcol, int* Brow, int** Ccol, int* Crow, int *Csize);


#endif


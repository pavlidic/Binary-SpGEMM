#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>

int rndCSR(int* row, int* col, int n, int m, float sparcity);

void full2CSR(const bool *C, const int n, const int m, int *Ccol, int *Crow);

bool full2CSR_re(const bool *C, const int n, const int m, int **Ccol, int *Crow, int *Csize);

void full2CSR_T(const bool *C, const int n, const int m, int *Ccol, int *Crow);

void printCSR(int* row, int* col, int n, int m, int b);


// function to swap elements
void swap(int *a, int *b);

// function to find the partition position
int partition(int array[], int low, int high);

void quickSort(int array[], int low, int high);



#endif


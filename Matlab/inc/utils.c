#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include "utils.h"

int rndCSR(int* row, int* col, int n, int m, float sparcity){
    float chance;
    int sum=0;

    row[0]=0;
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            chance = (float)rand()/RAND_MAX;
            if(chance<sparcity){
                col[sum]=j;
                sum++;
            }
        }
        row[i+1]=sum;
    }

    return row[n];
}

void full2CSR(const bool *C, const int n, const int m, int *Ccol, int *Crow){
    int nnz=0;
    for(int i=0; i<n; i++){
        Crow[i]=nnz;
        for(int j=0; j<m; j++){
            if(C[i*m + j]){
                Ccol[nnz++]=j;
            }
        }        
    }
    Crow[n]=nnz;
}

bool full2CSR_re(const bool *C, const int n, const int m, int **Ccol, int *Crow, int *Csize){
    int nnz=0;
    for(int i=0; i<n; i++){
        Crow[i]=nnz;
        for(int j=0; j<m; j++){
            if(C[i*m + j]){
                (*Ccol)[nnz++]=j;
                if(nnz==*Csize){
                    *Csize *= 2;
                    *Ccol = realloc(*Ccol,*Csize * sizeof(int));
                }
            }

        }  

    }
    Crow[n]=nnz;

    return nnz;
}

void full2CSR_T(const bool *C, const int n, const int m, int *Ccol, int *Crow){
    int nnz=0;
    for(int i=0; i<n; i++){
        Crow[i]=nnz;
        for(int j=0; j<m; j++){
            if(C[j*m + i]){
                Ccol[nnz++]=j;
            }
        }        
    }
    Crow[n]=nnz;
}

void printCSR(int* row, int* col, int n, int m, int b){
    int* temp_row = malloc((n+1)*sizeof(int));
    //fix row collumn if biased
    int bias = row[0];
    for(int i=0; i<n+1; i++){
        temp_row[i]=row[i]-bias;
    }

    int row_nnz_count = 0;
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            if(temp_row[i]+row_nnz_count < temp_row[i+1] && col[temp_row[i]+row_nnz_count] == j){
                printf("1 ");
                row_nnz_count++;
            }else{
                printf("0 ");
            }
            if((j+1)%b==0 && j+1<m){
                printf("| ");
            }
        }
        row_nnz_count=0;
        printf("\n");
        if((i+1)%b==0 && i+1<n ){
            for(int j=0; j<m + m/b - 1; j++){
                printf("--");
            }
            printf("\n");
        }
    }
    printf("\n");
}

// function to swap elements
void swap(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

// function to find the partition position
int partition(int array[], int low, int high) {
  
  // select the rightmost element as pivot
  int pivot = array[high];
  
  // pointer for greater element
  int i = (low - 1);

  // traverse each element of the array
  // compare them with the pivot
  for (int j = low; j < high; j++) {
    if (array[j] <= pivot) {
        
      // if element smaller than pivot is found
      // swap it with the greater element pointed by i
      i++;
      
      // swap element at i with element at j
      swap(&array[i], &array[j]);
    }
  }

  // swap the pivot element with the greater element at i
  swap(&array[i + 1], &array[high]);
  
  // return the partition point
  return (i + 1);
}

void quickSort(int array[], int low, int high) {
  if (low < high) {
    
    // find the pivot element such that
    // elements smaller than pivot are on left of pivot
    // elements greater than pivot are on right of pivot
    int pi = partition(array, low, high);
    
    // recursive call on the left of pivot
    quickSort(array, low, pi - 1);
    
    // recursive call on the right of pivot
    quickSort(array, pi + 1, high);
  }
}


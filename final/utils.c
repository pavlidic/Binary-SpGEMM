#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include "coo2csc.h"
#include "mmio.h"
#include "utils.h"




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

void readCOO(const char *mat, uint32_t **row, uint32_t **col,uint32_t *M, uint32_t *N, uint32_t *nnz){
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t *I, *J;
    uint32_t isOneBased=0;

    if ((f = fopen(mat, "r")) == NULL) 
        exit(1);
    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nnz)) !=0)
        exit(1);

    I = malloc(*nnz * sizeof(uint32_t));
    J = malloc(*nnz * sizeof(uint32_t));

    for (int i=0; i<*nnz; i++)
    {
        fscanf(f, "%u %u\n", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }
    if (f !=stdin) fclose(f);

    *row  = malloc((*M + 1) * sizeof(uint32_t));   
    *col  = malloc(*nnz  * sizeof(uint32_t));

    coo2csc(*col, *row, I,J,*nnz,*M,isOneBased);

    free(I);
    free(J);
}


struct timespec diff(struct timespec start, struct timespec end)
{
        struct timespec temp;
        if ((end.tv_nsec - start.tv_nsec) < 0) 
        {
                temp.tv_sec = end.tv_sec - start.tv_sec - 1;
                temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
        } 
        else 
        {
                temp.tv_sec = end.tv_sec - start.tv_sec;
                temp.tv_nsec = end.tv_nsec - start.tv_nsec;
        }
        return temp;
}

double timeConv(struct timespec result){
        return (double)(result.tv_sec+(double)result.tv_nsec/1000000000);
}

double tictoc(int mode) {
    static struct timespec t_start, t_end;
    
    if (mode==0)
        clock_gettime(CLOCK_MONOTONIC,&t_start);
    else {
        clock_gettime(CLOCK_MONOTONIC,&t_end);
        return timeConv(diff(t_start,t_end));
    }
}


// function to swap elements
void swap(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

void swapD(double *a, double *b) {
  double t = *a;
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

int partitionD(double array[], int low, int high) {
  
  // select the rightmost element as pivot
  double pivot = array[high];
  
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
      swapD(&array[i], &array[j]);
    }
  }

  // swap the pivot element with the greater element at i
  swapD(&array[i + 1], &array[high]);
  
  // return the partition point
  return (i + 1);
}

void quickSortD(double array[], int low, int high) {
  if (low < high) {
    
    // find the pivot element such that
    // elements smaller than pivot are on left of pivot
    // elements greater than pivot are on right of pivot
    int pi = partitionD(array, low, high);
    
    // recursive call on the left of pivot
    quickSortD(array, low, pi - 1);
    
    // recursive call on the right of pivot
    quickSortD(array, pi + 1, high);
  }
}






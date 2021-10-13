#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include "coo2csc.h"
#include "mmio.h"
#include "utils.h"


int rndCSR(int** row, int** col, int n, int m, int nnz_per_row){
    float chance;
    int sum=0;

    float sparcity = ((float)nnz_per_row)/n;

    *row = malloc((n+1)*sizeof(int));
    *col = malloc((2*nnz_per_row*n)*sizeof(int));

    (*row)[0]=0;
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            chance = (float)rand()/RAND_MAX;
            //printf("ch:%f, sp:%f\n",chance,sparcity);
            if(chance<sparcity){
                (*col)[sum]=j;
                sum++;
            }
        }
        (*row)[i+1]=sum;
    }
    //printf("sum:%d\n",sum);

    *col = realloc(*col,sum*sizeof(int));

    return (*row)[n];
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

//saves a CSR matrix to a blocked version of it
void csr2bcsr(int*  csr_row, int*  csr_col, int n, int m, int nnz, int b, 
              CSRbCSR* A, int print){
    
    A->b = b;
    A->n = n;
    A->m = m;
    A->nnz = nnz;

    

    int col;
    int xB, yB;
    int num_of_blocks = (n/b)*(n/b), blocks_per_blockrow = n/b;
    int *Bmat_nz = malloc(num_of_blocks*sizeof(int));
    for(int i=0; i<num_of_blocks; i++){
        Bmat_nz[i]=0;
    }

    A->Bcol = malloc(num_of_blocks*sizeof(int));
    A->Brow = malloc((blocks_per_blockrow+1)*sizeof(int));

    for(int i=0; i<n; i++){
        for(int j=csr_row[i]; j<csr_row[i+1]; j++){

            col = csr_col[j];
            xB = col/b;
            yB =   i/b;
            Bmat_nz[yB*blocks_per_blockrow + xB]++;

        }
    }

    int *Bmat_nz_sum = malloc(num_of_blocks*sizeof(int));
    Bmat_nz_sum[0] = Bmat_nz[0];
    for(int i=1; i<num_of_blocks; i++){
        Bmat_nz_sum[i] = Bmat_nz_sum[i-1] + Bmat_nz[i];
    }

    
    //create High-Level Block csr
    int nnzb_count=0; //count for blocks with at least 1 nz
    int *block_to_nz_block = malloc(num_of_blocks*sizeof(int)); //maps blocks to row-major blocks with at least 1 nz
                                                                //eg. block 5(4+1)-> nz block 3(2+1) if only 2 out of 4 previous blocks had a nz
    for(int i=0; i<n/b; i++){
        for(int j=0; j<n/b; j++){
            block_to_nz_block[i*blocks_per_blockrow + j] = -1;
        }
    }
    /* for(int i=0; i<n/b; i++){
        for(int j=0; j<n/b; j++){
            Bcsr_col[i*blocks_per_blockrow + j] = -1;
        }
    } */

    A->Brow[0] = 0;
    for(int i=0; i<n/b; i++){
        for(int j=0; j<n/b; j++){
            if(Bmat_nz[i*blocks_per_blockrow + j]!=0 ){//&& Bcsr_col[nnzb_count]==0){//if block has nonzero and hasnt been added yet -> add it
                block_to_nz_block[i*blocks_per_blockrow + j] = nnzb_count; 
                A->Bcol[nnzb_count++] = j;
            }
        }
        A->Brow[i+1] = nnzb_count;
    }
    A->nnzb = nnzb_count;

    
    //create Low-Level csr
    A->LLcol = malloc(nnz*sizeof(int));
    A->LLrow = malloc((b*nnzb_count + 1)*sizeof(int));
    A->LLrow[0]=0;
    

    int bNum; //block number

    int *block_row_nnz_count = malloc(n*sizeof(int)); //counts nnz in each row of a block //n = n/b(number of blocks in a blockrow) * b(rows in a block)
    for(int j=0; j<n; j++){ block_row_nnz_count[j]=0;}

    int row_in_block, block_in_blockrow, current_block;

    for(int i=0; i<n; i++){

        for(int j=csr_row[i]; j<csr_row[i+1]; j++){

            col = csr_col[j];
            xB = col/b;
            yB =   i/b;
            row_in_block = i%b;

            bNum = yB*blocks_per_blockrow + xB; 

            A->LLcol[ Bmat_nz_sum[bNum]-Bmat_nz[bNum] ] = col%b; //fill LLcol from the start of each block
            Bmat_nz[bNum]--; //go to next element of each block in LLcol
            block_row_nnz_count[xB*b + row_in_block]++;
        }

        if((i+1)%b==0){ 
            for(int j=0; j<n; j++){
                yB =                i/b;
                block_in_blockrow = j/b;
                row_in_block =      j%b;
                current_block = yB*blocks_per_blockrow+block_in_blockrow;
                if(block_to_nz_block[current_block]!=-1){ //if at a block that has a nz (maps to a nz block)
                    A->LLrow[ block_to_nz_block[current_block]*b + row_in_block + 1 ] = block_row_nnz_count[j] + A->LLrow[ block_to_nz_block[current_block]*b + row_in_block ];
                }
                block_row_nnz_count[j]=0;
            }
        }
    }

    free(Bmat_nz);
    free(block_row_nnz_count);
    free(block_to_nz_block);
    free(Bmat_nz_sum);


    //print
    if(print){
        printf("Bcol = ");
        for(int i=0; i<nnzb_count; i++){
            printf("%d, ",A->Bcol[i]);
        }
        printf("\nBrow = ");
        for(int i=0; i<n/b+1; i++){
            printf("%d, ",A->Brow[i]);
        }
        printf("\nLLcol = ");
        for(int i=0; i<nnz; i++){
            printf("%d, ",A->LLcol[i]);
        }
        printf("\nLLrow = ");
        for(int i=0; i<b*nnzb_count+1; i++){
            printf("%d, ",A->LLrow[i]);
        }
        printf("\n");
    }
    

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

//TODO:test
/* void quickSort(int * tab, int l, int r)
{
   int q;
   while(l < r)
   {
      q = partition(tab, l, r);
      if(q - l < r - q) //recurse into the smaller half
      {
         quickSort(tab, l, q - 1);
         l = q + 1;
      } else
      {
         quickSort(tab, q + 1, r);
         r = q - 1;
      }
   }
} */


//merges 2 arrays into a third, no doubles, returns n3
int merge(int* arr1, int* arr2, int* arr3, int n1, int n2){
    int i = 0, j = 0, k = 0;
    int n3=0;

    while(i<n1 && j <n2)
    {
        if(arr1[i] < arr2[j])
            arr3[k++] = arr1[i++];
        else if(arr1[i] > arr2[j])
            arr3[k++] = arr2[j++];
        else{
            arr3[k++] = arr2[j++];
            i++;
        }
        n3++;
    }
 
    // Store remaining elements of first array
    while (i < n1){
        arr3[k++] = arr1[i++];
        n3++;
    }
        
 
    // Store remaining elements of second array
    while (j < n2){
        arr3[k++] = arr2[j++];
        n3++;
    }
    
    return n3;
}

// ORs 2 sparse matrices(CSR)
void SpM_OR(int* Acol, int* Arow, int n,
            int* Bcol, int* Brow,
            int** Ccol, int* Crow, int *Csize)
{
    if(*Csize < Arow[n]+Brow[n]){ //check if enough space in Ccol
        *Csize = Arow[n]+Brow[n];
        *Ccol = realloc(*Ccol, *Csize*sizeof(int));
    }

    int cumsum = 0;
    Crow[0] = 0;
    for(int i=0; i<n; i++){
        cumsum += merge(&Acol[Arow[i]], &Bcol[Brow[i]], &(*Ccol)[cumsum], Arow[i+1]-Arow[i], Brow[i+1]-Brow[i]);
        Crow[i+1] = cumsum;
    }

}


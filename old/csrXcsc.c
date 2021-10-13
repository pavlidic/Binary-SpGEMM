#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include "timediff.h"
#include "coo2csc.h"
#include "mmio.h"

//TODO: ORing of submatrices with previous ones maybe happens while calculating the BMM


typedef struct CSRbCSR{
    int  n,m;
    int* LLcol;
    int* LLrow;
    int* Bcol;
    int* Brow;
    int  nnz;
    int  nnzb;
    int  b;
}CSRbCSR;

// creates a random sparse matrix in CSR form, returns NNZ
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


//prints a binary CSR matrix, in b-length blocks. if no blocks, set b=n
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

//saves a CSR matrix to a blocked version of it
void csr2bcsr(int*  csr_row, int*  csr_col, int n, int nnz, int b, 
              CSRbCSR* A, int print){
    
    A->b = b;
    A->n = n; //TODO:add m
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

//Boolean sparse Matrix Multiplication
void BSpMM(int* Acol, int* Arow, int An, int Annz,
           int* Bcol, int* Brow, int Bn, int Bnnz,
           int print)
{
    int Crow[An+1];
    int Ccol[1000];//Cnnz
    int Crow_sum = 0; // to immediately produce C in csr 

    for(int i = 0; i<An; i++){

        Crow[i]=Crow_sum;

        for(int j = 0; j<Bn; j++){

            int t1 = Arow[i];
            int t2 = Bcol[j];
            int n1 = Arow[i+1];
            int n2 = Bcol[j+1];

            while(t1<n1 && t2<n2){
                if(Acol[t1]==Brow[t2]){
                    Ccol[Crow_sum]=j;
                    Crow_sum++;
                    break;
                }else if(Acol[t1]<Brow[t2]){
                    t1++;
                }else{
                    t2++;
                }
            }
        }        
    }
    Crow[An] = Crow_sum;

    if(print==1){
        printf("\n\nCcol= ");
        for(int i=0; i<Crow_sum; i++){
            printf("%d ",Ccol[i]);
        }
        printf("\nCrow= ");
        for(int i=0; i<An+1; i++){
            printf("%d ",Crow[i]);
        }
    }else if(print==2){
        printf("C=\n");
        printCSR(Crow,Ccol,An,An,An);
    }

}

//Boolean sparse Matrix Multiplication
bool BSpMM2(int* Acol, int* Arow, int An, // returns 1 if nzb, 0 if zb
            int* Bcol, int* Brow, int Bn,
            int* Ccol, int* Crow, int Cn,
            int print)
{

    int Crow_sum = 0; // to immediately produce C in csr 
    Crow = malloc((An+1)*sizeof(int));
    int Ccol_size = An; //initial size
    Ccol = malloc(Ccol_size*sizeof(int));

    for(int i = 0; i<An; i++){

        Crow[i]=Crow_sum;

        for(int j = 0; j<Bn; j++){

            int t1 = Arow[i];
            int t2 = Bcol[j];
            int n1 = Arow[i+1];
            int n2 = Bcol[j+1];

            while(t1<n1 && t2<n2){ // for C(i,j) check A(i,k)*B(k,j) for each k
                if(Acol[t1]==Brow[t2]){
                    Ccol[Crow_sum]=j;
                    Crow_sum++;
                    if(Crow_sum==Ccol_size){ //check if you need more space
                        Ccol_size *= 2;
                        Ccol = realloc(Ccol, Ccol_size*sizeof(int));
                    }
                    break;
                }else if(Acol[t1]<Brow[t2]){
                    t1++;
                }else{
                    t2++;
                }
            }

        }        
    }
    Crow[An] = Crow_sum;

    //print
    if(print==1){
        printf("\n\nCcol= ");
        for(int i=0; i<Crow_sum; i++){
            printf("%d ",Ccol[i]);
        }
        printf("\nCrow= ");
        for(int i=0; i<An+1; i++){
            printf("%d ",Crow[i]);
        }
    }else if(print==2){
        printf("C=\n");
        printCSR(Crow,Ccol,An,An,An);
    }


    if(Crow_sum){
        Ccol = realloc(Ccol, Crow_sum*sizeof(int)); //use just enough space
        return 1;
    }else{
        free(Ccol);
        free(Crow);
        return 0;
    }
}

//Boolean sparse Matrix Multiplication
bool BSpMM3(int*  Acol, int* Arow, int An,                  // returns 1 if nzb, 0 if zb
            int*  Bcol, int* Brow, int Bn,
            int** Ccol, int* Crow, int Cn, int *Ccol_size,
            int print)
{

    int Crow_sum = 0; // to immediately produce C in csr 
    /* Crow = malloc((An+1)*sizeof(int));
    int Ccol_size = An; //initial size
    Ccol = malloc(Ccol_size*sizeof(int)); */

    uint32_t i,j,t1,t2,n1,n2;

    for(i = 0; i<An; i++){

        Crow[i]=Crow_sum;

        for(j = 0; j<Bn; j++){

            t1 = Arow[i];
            t2 = Bcol[j];
            n1 = Arow[i+1];
            n2 = Bcol[j+1];

            while(t1<n1 && t2<n2){ // for C(i,j) check A(i,k)*B(k,j) for each k
                if(Acol[t1]==Brow[t2]){
                    (*Ccol)[Crow_sum]=j;
                    Crow_sum++;
                    if(Crow_sum == *Ccol_size){ //check if you need more space
                        *Ccol_size *= 2;
                        *Ccol = realloc(*Ccol, *Ccol_size*sizeof(int));
                    }
                    break;
                }else if(Acol[t1]<Brow[t2]){
                    t1++;
                }else{
                    t2++;
                }
            }

        }        
    }
    Crow[An] = Crow_sum;

    //print
    if(print==1){
        printf("\n\nCcol= ");
        for(int i=0; i<Crow_sum; i++){
            printf("%d ",*Ccol[i]);
        }
        printf("\nCrow= ");
        for(int i=0; i<An+1; i++){
            printf("%d ",Crow[i]);
        }
    }else if(print==2){
        printf("C=\n");
        printCSR(Crow,*Ccol,An,An,An);
    }


    return Crow_sum;
}

bool BSpMM4(int*  Acol, int* Arow, int An,                  // returns 1 if nzb, 0 if zb
            int*  Bcol, int* Brow, int Bn,                  // both CSR
            int** Ccol, int* Crow, int Cn, int *Ccol_size,
            int print)
{
    bool *C = calloc(An*Bn,sizeof(bool));
    bool ret = 0;
    for(int i=0; i<An; i++){
        for(int jj=Arow[i]; jj<Arow[i+1]; jj++){
            int j = Acol[jj];

            for(int k=Brow[j]; k<Brow[j+1]; k++){
                C[i*An + Bcol[k]] = 1;
            }

        }
    }

    ret = full2CSR_re(C,An,Bn,Ccol,Crow,Ccol_size);
    free(C);
    return ret;

}

bool BSpMM_masked(int*  Acol, int* Arow, int An,                  // returns 1 if nzb, 0 if zb
                  int*  Bcol, int* Brow, int Bn,
                  int*  Fcol, int* Frow, int Fn,
                  int** Ccol, int* Crow, int Cn, int *Ccol_size,
                  int print)
{

    int Crow_sum = 0; // to immediately produce C in csr 

    for(int i = 0; i<Fn; i++){

        Crow[i]=Crow_sum;

        for(int k = Frow[i]; k<Frow[i+1]; k++){

            int j = Fcol[k];

            int t1 = Arow[i];
            int t2 = Bcol[j];
            int n1 = Arow[i+1];
            int n2 = Bcol[j+1];

            while(t1<n1 && t2<n2){ // for C(i,j) check A(i,k)*B(k,j) for each k
                if(Acol[t1]==Brow[t2]){
                    (*Ccol)[Crow_sum]=j;
                    Crow_sum++;
                    if(Crow_sum == *Ccol_size){ //check if you need more space
                        *Ccol_size *= 2;
                        *Ccol = realloc(*Ccol, *Ccol_size*sizeof(int));
                    }
                    break;
                }else if(Acol[t1]<Brow[t2]){
                    t1++;
                }else{
                    t2++;
                }
            }

        }        
    }
    Crow[An] = Crow_sum;

    //print
    if(print==1){
        printf("\n\nCcol= ");
        for(int i=0; i<Crow_sum; i++){
            printf("%d ",*Ccol[i]);
        }
        printf("\nCrow= ");
        for(int i=0; i<An+1; i++){
            printf("%d ",Crow[i]);
        }
    }else if(print==2){
        printf("C=\n");
        printCSR(Crow,*Ccol,An,An,An);
    }


    if(Crow_sum){
        return 1;
    }else{
        return 0;
    }
}


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

void SpBMM(int* Acol, int* Arow, int* ABcol, int* ABrow,int An, int Annz, int Annzb, //TODO: currently only square MM, change to whatever
           int* Bcol, int* Brow, int* BBcol, int* BBrow,int Bn, int Bnnz, int Bnnzb,
           int* Ccol, int* Crow, int* CBcol, int* CBrow,int Cn, int Cnnz, int Cnnzb,
           int b, int print){
    
    int bpr = An/b; //blocks per row
    for(int i=0; i<bpr; i++){
        for(int j=ABrow[i]; j<ABrow[i+1]; j++){

        }
    }

    int* Crow_tmp, Ccol_tmp;
    for(int i=0; i<bpr; i++){       // i = C's blockRow ptr
        for(int j=0; j<bpr; j++){   // j = C's blockCol ptr

            int t1 = ABrow[i];
            int t2 = BBcol[j];
            int n1 = ABrow[i+1];
            int n2 = BBcol[j+1];

            while(t1<n1 && t2<n2){
                if(ABcol[t1]==BBrow[t2]){//do MM

                    /* int xA = b*ABcol[t1];
                    int xB = b*BBrow[t2];

                    int yA = b*i;
                    int yB = b*j; */

                    int Aptr = b*t1; // points to the start of the block in LLrow
                    int Bptr = b*t2;
                    
                    int Ablock_nnz = Arow[Aptr+b]-Arow[Aptr];           //TODO: get directly from bcsr creation, also not needed
                    int Bblock_nnz = Bcol[Bptr+b]-Bcol[Bptr];
                    

                    printf("A(%d,%d)*B(%d,%d)\n",i,ABcol[t1],BBrow[t2],j);

                    printCSR(&Arow[Aptr],&Acol[ Arow[Aptr] ],b,b,b);    //FIXME:these dont really work cuase 1 is csr other is csc
                    printCSR(&Bcol[Bptr],&Brow[ Bcol[Bptr] ],b,b,b);    //<- should be printCSC(could convert between CSR and CSC and print with the same function)

                    BSpMM(Acol, &Arow[Aptr], b, Ablock_nnz,             //TODO: maybe change it so you send a ptr to Acol where it starts (need to change how LLrow works)
                          &Bcol[Bptr], Brow, b, Bblock_nnz,
                          print);

                    t1++;
                    t2++;

                }else if(ABcol[t1]<BBrow[t2]){
                    t1++;
                }else{
                    t2++;
                }
            }

        }
    }
}

void SpBMM2(int*  Acol, int* Arow, int* ABcol, int* ABrow,int An, int Annz, int Annzb, //TODO: currently only square MM, change to whatever
            int*  Bcol, int* Brow, int* BBcol, int* BBrow,int Bn, int Bnnz, int Bnnzb,
            int** Ccol, int* Crow, int* CBcol, int* CBrow,int Cn, int *Cnnz,
            int b, int print){
    
    int bpr = An/b; //blocks per row
    int Abpc = An/b; //blocks per collumn

    int init_size = Abpc;

    
    int NB = Abpc*Abpc; //number of blocks
    int  *Ccol_size1 = malloc(NB*sizeof(int)); //contains the sizes for Ccol of each block //TODO: free
    int  *Ccol_size2 = malloc(NB*sizeof(int));
    int  *Ccol_size3 = malloc(NB*sizeof(int));

    int **Ccol_dist1 = malloc(NB*sizeof(int*));//distributed
    int **Ccol_dist2 = malloc(NB*sizeof(int*));// 2 |= 1
    int **Ccol_dist3 = malloc(NB*sizeof(int*));

    int **Crow_dist1 = malloc(NB*sizeof(int*));//distributed
    int **Crow_dist2 = malloc(NB*sizeof(int*));// 2 |= 1
    int **Crow_dist3 = malloc(NB*sizeof(int*));

    int* Crow_0 = malloc((b+1)*sizeof(int));
    for(int i=0; i<(b+1); i++){ Crow_0[i] = 0; }

    bool *is_nzb = malloc(NB*sizeof(bool));
    for(int i=0; i<NB; i++){ is_nzb[i] = 0; }

    for(int i=0; i<NB; i++){
        Ccol_size1[i] = init_size;
        Ccol_size2[i] = init_size;
        Ccol_size3[i] = init_size;
    }

    //foreach block MM
    for(int i=0; i<bpr; i++){       // i = C's blockRow ptr
        for(int j=0; j<bpr; j++){   // j = C's blockCol ptr
            bool first_time = 1;
            int blk = i*bpr + j;

            bool cnt = 0;

            int t1 = ABrow[i];
            int t2 = BBcol[j];
            int n1 = ABrow[i+1];
            int n2 = BBcol[j+1];

            while(t1<n1 && t2<n2){
                if(ABcol[t1]==BBrow[t2]){//do MM

                    int Aptr = b*t1; // points to the start of the block in LLrow
                    int Bptr = b*t2;

                    // C matrix allocation
                    if(first_time){
                        first_time = 0;
                        Ccol_dist1[blk] = malloc(init_size*sizeof(int));
                        Ccol_dist2[blk] = malloc(init_size*sizeof(int));
                        Ccol_dist3[blk] = malloc(init_size*sizeof(int));
                        Crow_dist1[blk] = malloc((b+1)*sizeof(int));
                        Crow_dist2[blk] = malloc((b+1)*sizeof(int));
                        Crow_dist3[blk] = malloc((b+1)*sizeof(int));
                        memcpy(Crow_dist2[blk], Crow_0, (b+1)*sizeof(int)); // initially zero matrix
                    }

                    //printf("A(%d,%d)*B(%d,%d)\n",i,ABcol[t1],BBrow[t2],j);
                    //printCSR(&Arow[Aptr],&Acol[ Arow[Aptr] ],b,b,b);    //FIXME:these dont really work cuase 1 is csr other is csc
                    //printCSR(&Bcol[Bptr],&Brow[ Bcol[Bptr] ],b,b,b);    //<- should be printCSC(could convert between CSR and CSC and print with the same function)

                    BSpMM3(Acol, &Arow[Aptr], b,             //TODO: maybe change it so you send a ptr to Acol where it starts (need to change how LLrow works)
                           &Bcol[Bptr], Brow, b,
                           &Ccol_dist1[blk], Crow_dist1[blk], b, &Ccol_size1[blk],
                           print);
                    
                    //printCSR(Crow_dist1[blk],Ccol_dist1[blk],b,b,b);
                    
                    //printf("ORed=\n");
                    if(cnt==0){ //alternate between Ccol_dist 2 and 3 to save a memcpy
                        SpM_OR(  Ccol_dist1[blk], Crow_dist1[blk], b,
                                 Ccol_dist2[blk], Crow_dist2[blk],
                                &Ccol_dist3[blk], Crow_dist3[blk], &Ccol_size3[blk]);
                        cnt=1;
                        //printCSR(Crow_dist3[blk],Ccol_dist3[blk],b,b,b);
                    }else{
                        SpM_OR(  Ccol_dist1[blk], Crow_dist1[blk], b,
                                 Ccol_dist3[blk], Crow_dist3[blk],
                                &Ccol_dist2[blk], Crow_dist2[blk], &Ccol_size2[blk]);
                        cnt=0;
                        //printCSR(Crow_dist2[blk],Ccol_dist2[blk],b,b,b);
                    }

                    t1++;
                    t2++;
                }else if(ABcol[t1]<BBrow[t2]){
                    t1++;
                }else{
                    t2++;
                }
            }



            if(!first_time){ // change pointers so that final matrix is in C_dist1, free space
                int *temp_col_ptr;
                int *temp_row_ptr;

                if(cnt==0){
                    temp_col_ptr = Ccol_dist1[blk];
                    Ccol_dist1[blk] = Ccol_dist2[blk];
                    Ccol_dist2[blk] = temp_col_ptr;

                    temp_row_ptr = Crow_dist1[blk];
                    Crow_dist1[blk] = Crow_dist2[blk];
                    Crow_dist2[blk] = temp_row_ptr;
                }else{
                    temp_col_ptr = Ccol_dist1[blk];
                    Ccol_dist1[blk] = Ccol_dist3[blk];
                    Ccol_dist3[blk] = temp_col_ptr;

                    temp_row_ptr = Crow_dist1[blk];
                    Crow_dist1[blk] = Crow_dist3[blk];
                    Crow_dist3[blk] = temp_row_ptr;
                }

                free(Ccol_dist2[blk]);
                free(Ccol_dist3[blk]);

                free(Crow_dist2[blk]);
                free(Crow_dist3[blk]);

                if(Crow_dist1[blk][b]) is_nzb[blk] = 1;
            }

        }
    }
    printf("blocked calc done\n");


    *Cnnz = 0;
    for(int i=0; i<Abpc; i++){
        for(int j=0; j<Abpc; j++){
            int blk = i*Abpc + j;
            if(is_nzb[blk]){
                *Cnnz += Crow_dist1[blk][b];
            }
        }
    }

    *Ccol = malloc(*Cnnz*sizeof(int));

    int nzb_cnt = 0;
    //reconstruct C
    Crow[0] = 0;
    for(int iB=0; iB<Abpc; iB++){
        for(int i=0; i<b; i++){
            for(int jB=0; jB<Abpc; jB++){
                int blk = iB*Abpc + jB;
                //printCSR(Crow_dist1[blk],Ccol_dist1[blk],b,b,b);
                if(is_nzb[blk]){
                    for(int j=Crow_dist1[blk][i]; j<Crow_dist1[blk][i+1]; j++){
                        (*Ccol)[nzb_cnt] = jB*b + Ccol_dist1[blk][j];
                        nzb_cnt++;
                    }

                    if( i==(b-1) ){//if last row in that block, free block
                        free(Crow_dist1[blk]);
                        free(Ccol_dist1[blk]);
                    }
                }

                
            }

            Crow[iB*b + i + 1] = nzb_cnt;
        }
    }
    free(Ccol_size1);
    free(Ccol_size2);
    free(Ccol_size3);
    free(Ccol_dist1);
    free(Ccol_dist2);
    free(Ccol_dist3);
    free(Crow_dist1);
    free(Crow_dist2);
    free(Crow_dist3);
    free(is_nzb);
    free(Crow_0);

    //printCSR(Crow,*Ccol,An,An,b);

}

void SpBMM3(int*  Acol, int* Arow, int* ABcol, int* ABrow,int An, int Annz, int Annzb, //TODO: currently only square MM, change to whatever
            int*  Bcol, int* Brow, int* BBcol, int* BBrow,int Bn, int Bnnz, int Bnnzb,
            int** Ccol, int* Crow, int* CBcol, int* CBrow,int Cn, int *Cnnz,
            int b, int print){
    
    int bpr = An/b; //blocks per row
    int Abpc = An/b; //blocks per collumn

    int init_size = Abpc;

    
    int NB = Abpc*Abpc; //number of blocks
    int  *Ccol_size1 = malloc(NB*sizeof(int)); //contains the sizes for Ccol of each block //TODO: free
    int  *Ccol_size2 = malloc(NB*sizeof(int));
    int  *Ccol_size3 = malloc(NB*sizeof(int));

    int **Ccol_dist1 = malloc(NB*sizeof(int*));//distributed
    int **Ccol_dist2 = malloc(NB*sizeof(int*));// 2 |= 1
    int **Ccol_dist3 = malloc(NB*sizeof(int*));

    int **Crow_dist1 = malloc(NB*sizeof(int*));//distributed
    int **Crow_dist2 = malloc(NB*sizeof(int*));// 2 |= 1
    int **Crow_dist3 = malloc(NB*sizeof(int*));

    int* Crow_0 = malloc((b+1)*sizeof(int));
    for(int i=0; i<(b+1); i++){ Crow_0[i] = 0; }

    bool *is_nzb = malloc(NB*sizeof(bool));
    for(int i=0; i<NB; i++){ is_nzb[i] = 0; }

    for(int i=0; i<NB; i++){
        Ccol_size1[i] = init_size;
        Ccol_size2[i] = init_size;
        Ccol_size3[i] = init_size;
    }

    //foreach block MM
    for(int i=0; i<bpr; i++){       // i = C's blockRow ptr
        for(int j=0; j<bpr; j++){   // j = C's blockCol ptr
            bool first_time = 1;
            int blk = i*bpr + j;

            bool cnt = 0;

            int t1 = ABrow[i];
            int t2 = BBrow[j];
            int n1 = ABrow[i+1];
            int n2 = BBrow[j+1];

            while(t1<n1 && t2<n2){
                if(ABcol[t1]==BBcol[t2]){//do MM

                    int Aptr = b*t1; // points to the start of the block in LLrow
                    int Bptr = b*t2;

                    // C matrix allocation
                    if(first_time){
                        first_time = 0;
                        Ccol_dist1[blk] = malloc(init_size*sizeof(int));
                        Ccol_dist2[blk] = malloc(init_size*sizeof(int));
                        Ccol_dist3[blk] = malloc(init_size*sizeof(int));
                        Crow_dist1[blk] = malloc((b+1)*sizeof(int));
                        Crow_dist2[blk] = malloc((b+1)*sizeof(int));
                        Crow_dist3[blk] = malloc((b+1)*sizeof(int));
                        memcpy(Crow_dist2[blk], Crow_0, (b+1)*sizeof(int)); // initially zero matrix
                    }

                    //printf("A(%d,%d)*B(%d,%d)\n",i,ABcol[t1],BBrow[t2],j);
                    //printCSR(&Arow[Aptr],&Acol[ Arow[Aptr] ],b,b,b);    //FIXME:these dont really work cuase 1 is csr other is csc
                    //printCSR(&Bcol[Bptr],&Brow[ Bcol[Bptr] ],b,b,b);    //<- should be printCSC(could convert between CSR and CSC and print with the same function)

                    BSpMM4(Acol, &Arow[Aptr], b,             //TODO: maybe change it so you send a ptr to Acol where it starts (need to change how LLrow works)
                           Bcol, &Brow[Bptr], b,
                           &Ccol_dist1[blk], Crow_dist1[blk], b, &Ccol_size1[blk],
                           print);
                    
                    //printCSR(Crow_dist1[blk],Ccol_dist1[blk],b,b,b);
                    
                    //printf("ORed=\n");
                    if(cnt==0){ //alternate between Ccol_dist 2 and 3 to save a memcpy
                        SpM_OR(  Ccol_dist1[blk], Crow_dist1[blk], b,
                                 Ccol_dist2[blk], Crow_dist2[blk],
                                &Ccol_dist3[blk], Crow_dist3[blk], &Ccol_size3[blk]);
                        cnt=1;
                        //printCSR(Crow_dist3[blk],Ccol_dist3[blk],b,b,b);
                    }else{
                        SpM_OR(  Ccol_dist1[blk], Crow_dist1[blk], b,
                                 Ccol_dist3[blk], Crow_dist3[blk],
                                &Ccol_dist2[blk], Crow_dist2[blk], &Ccol_size2[blk]);
                        cnt=0;
                        //printCSR(Crow_dist2[blk],Ccol_dist2[blk],b,b,b);
                    }

                    t1++;
                    t2++;
                }else if(ABcol[t1]<BBcol[t2]){
                    t1++;
                }else{
                    t2++;
                }
            }



            if(!first_time){ // change pointers so that final matrix is in C_dist1, free space
                int *temp_col_ptr;
                int *temp_row_ptr;

                if(cnt==0){
                    temp_col_ptr = Ccol_dist1[blk];
                    Ccol_dist1[blk] = Ccol_dist2[blk];
                    Ccol_dist2[blk] = temp_col_ptr;

                    temp_row_ptr = Crow_dist1[blk];
                    Crow_dist1[blk] = Crow_dist2[blk];
                    Crow_dist2[blk] = temp_row_ptr;
                }else{
                    temp_col_ptr = Ccol_dist1[blk];
                    Ccol_dist1[blk] = Ccol_dist3[blk];
                    Ccol_dist3[blk] = temp_col_ptr;

                    temp_row_ptr = Crow_dist1[blk];
                    Crow_dist1[blk] = Crow_dist3[blk];
                    Crow_dist3[blk] = temp_row_ptr;
                }

                free(Ccol_dist2[blk]);
                free(Ccol_dist3[blk]);

                free(Crow_dist2[blk]);
                free(Crow_dist3[blk]);

                if(Crow_dist1[blk][b]) is_nzb[blk] = 1;
            }

        }
    }
    printf("blocked calc done\n");


    *Cnnz = 0;
    for(int i=0; i<Abpc; i++){
        for(int j=0; j<Abpc; j++){
            int blk = i*Abpc + j;
            if(is_nzb[blk]){
                *Cnnz += Crow_dist1[blk][b];
            }
        }
    }

    *Ccol = malloc(*Cnnz*sizeof(int));

    int nzb_cnt = 0;
    //reconstruct C
    Crow[0] = 0;
    for(int iB=0; iB<Abpc; iB++){
        for(int i=0; i<b; i++){
            for(int jB=0; jB<Abpc; jB++){
                int blk = iB*Abpc + jB;
                //printCSR(Crow_dist1[blk],Ccol_dist1[blk],b,b,b);
                if(is_nzb[blk]){
                    for(int j=Crow_dist1[blk][i]; j<Crow_dist1[blk][i+1]; j++){
                        (*Ccol)[nzb_cnt] = jB*b + Ccol_dist1[blk][j];
                        nzb_cnt++;
                    }

                    if( i==(b-1) ){//if last row in that block, free block
                        free(Crow_dist1[blk]);
                        free(Ccol_dist1[blk]);
                    }
                }

                
            }

            Crow[iB*b + i + 1] = nzb_cnt;
        }
    }
    free(Ccol_size1);
    free(Ccol_size2);
    free(Ccol_size3);
    free(Ccol_dist1);
    free(Ccol_dist2);
    free(Ccol_dist3);
    free(Crow_dist1);
    free(Crow_dist2);
    free(Crow_dist3);
    free(is_nzb);
    free(Crow_0);

    //printCSR(Crow,*Ccol,An,An,b);

}


bool compare_CSR(int *Acol, int *Arow, int *Bcol, int *Brow, int n, bool print){
    for(int i=0; i<n+1; i++){
        if(print) printf("%d: Arow[%d]=%d, Brow[%d]=%d\n",i,i,Arow[i],i,Brow[i]);
        if(Arow[i]!=Brow[i]) {
            return 0;
        }
    }

    for(int i=0; i<Arow[n]; i++){
        if(print) printf("%d: Acol[%d]=%d, Bcol[%d]=%d\n",i,i,Acol[i],i,Bcol[i]);
        if(Acol[i]!=Bcol[i]) {
            return 0;
        }
    }

    return 1;
}

int test1(int argc, char const *argv[])
{
    srand(time(NULL));

    int An = 6, Am = 3;
    int Acol[8] = {0,1,2,0,2,0,1,0};
    int Arow[7] = {0,2,3,5,7,7,8};

    int Bn = 3, Bm = 8;
    int Brow[11] = {1,1,2,0,1,2,1,0,2,0,1};
    int Bcol[9] = {0,1,3,6,6,7,9,10,11};


    int Crow[1000];//An+1
    int Ccol[1000];//nnz
    int Crow_sum = 0; // to immediately produce C in csr 

    for(int i = 0; i<An; i++){

        Crow[i]=Crow_sum;

        for(int j = 0; j<Bm; j++){

            int t1 = Arow[i];
            int t2 = Bcol[j];
            int n1 = Arow[i+1];
            int n2 = Bcol[j+1];

            while(t1<n1 && t2<n2){
                if(Acol[t1]==Brow[t2]){
                    Ccol[Crow_sum]=j;
                    Crow_sum++;
                    break;
                }else if(Acol[t1]<Brow[t2]){
                    t1++;
                }else{
                    t2++;
                }
            }
        }        
    }
    Crow[An] = Crow_sum;

    printf("Ccol= ");
    for(int i=0; i<Crow_sum; i++){
        printf("%d ",Ccol[i]);
    }
    printf("\n\nCrow= ");
    for(int i=0; i<An+1; i++){
        printf("%d ",Crow[i]);
    }
     
    return 0;
}

// testing for when csr2bcsr didnt output to a CSRbCSR struct
/* int test2(int argc, char const *argv[])
{
    int col[] = {1,7,2,6,2,3,5,7,8,0,1,4,3};
    int row[] = {0,2,4,5,6,9,9,10,12,13};
    int n=9, b=3, nnz=13, nb = n*n/b/b;

    int Bcol[nb];
    int Brow[n/b+1];
    for(int i=0; i<nb; i++){
        Bcol[i]=0;
    }
    for(int i=0; i<n/b+1; i++){
        Brow[i]=0;
    }

    csr2bcsr(row,col,n,nnz,b,Brow,Bcol);

    
    return 0;
}
 */

void test3(){
    srand(time(NULL));

    int An=15, Am=15, Bn=15, Bm=15;
    int Arow[An+1],Brow[Bn+1];
    int Acol[1000], Bcol[1000];

    int Annz = rndCSR(Arow,Acol,An,Am,0.05);
    int Bnnz = rndCSR(Brow,Bcol,Bn,Bm,0.05);


    printf("\nAcol= ");
    for(int i=0; i<Annz; i++){
        printf("%d ",Acol[i]);
    }
    printf("\nArow= ");
    for(int i=0; i<An+1; i++){
        printf("%d ",Arow[i]);
    }

    printf("\n\nBcol= ");
    for(int i=0; i<Bnnz; i++){
        printf("%d ",Bcol[i]);
    }
    printf("\nBrow= ");
    for(int i=0; i<Bn+1; i++){
        printf("%d ",Brow[i]);
    }

    BSpMM(Acol,Arow,An,Annz,Brow,Bcol,Bn,Bnnz,1); //swap Brow/Bcol to make csr->csc transpose



}

void test4(){
    int col[] = {1,7,2,6,2,3,5,7,8,0,1,4,3};
    int row[] = {0,2,4,5,6,9,9,10,12,13};
    int n = 9;
    int nnz = 13;
    int b=3;

    //printCSR(row,col,n,n,b);

    CSRbCSR *A = malloc(sizeof(CSRbCSR));
    CSRbCSR *C = malloc(sizeof(CSRbCSR));

    csr2bcsr(row,col,n,nnz,b,A,1);

    SpBMM(A->LLcol, A->LLrow, A->Bcol, A->Brow, A->n, A->nnz, A->nnzb,
          A->LLrow, A->LLcol, A->Brow, A->Bcol, A->n, A->nnz, A->nnzb,
          C->LLcol, C->LLrow, C->Bcol, C->Brow, C->n, C->nnz, C->nnzb, b, 2);

}

void test5(){
    int arr1[] = {1,4,6,8,9,12};
    int arr2[] = {-2,2,3,6,8,9,10,11,12,13};
    int n1=6,n2=10;
    int* arr3 = malloc((n1+n2)*sizeof(int));

    int n3 = merge(arr1,arr2,arr3,n1,n2);

    for(int i=0; i<n3; i++){
        printf("%d ", arr3[i]);
    }
}

void test6(){
    int n=4;

    int Acol[] = {1,3,2,3,0,2,1};
    int Arow[] = {0,2,4,6,7};

    int Bcol[] = {0,1,0,0,1,2,3};
    int Brow[] = {0,2,3,5,7};

    int Csize = 7;
    int* Ccol = malloc(Csize*sizeof(int));
    int* Crow = malloc((n+1)*sizeof(int));

    SpM_OR(Acol,Arow,n,Bcol,Brow,&Ccol,Crow,&Csize);

    printCSR(Arow,Acol,n,n,n);
    printCSR(Brow,Bcol,n,n,n);
    printCSR(Crow,Ccol,n,n,n);
}

void test7(){
    int n=4;

    int Acol[] = {1,3,2,3,0,2,1};
    int Arow[] = {0,2,4,6,7};


    int* Ccol_0 = malloc(100*sizeof(int));
    int* Crow_0 = malloc((n+1)*sizeof(int));
    for(int i=0; i<(n+1); i++){ Crow_0[i] = 0; }

    int* Ccol = malloc(100*sizeof(int));
    int* Crow = malloc((n+1)*sizeof(int));
    int size = 100;
    SpM_OR(Acol,Arow,n,Ccol_0,Crow_0,&Ccol,Crow,&size);

    printCSR(Arow,Acol,n,n,n);
    printCSR(Crow_0,Ccol_0,n,n,n);
    printCSR(Crow,Ccol,n,n,n);

}

void test8(){
    int col[] = {1,7,2,6,2,3,5,7,8,0,1,4,3};
    int row[] = {0,2,4,5,6,9,9,10,12,13};
    int n = 9;
    int nnz = 13;
    int b=3;

    //printCSR(row,col,n,n,b);

    CSRbCSR *A = malloc(sizeof(CSRbCSR));
    CSRbCSR *C = malloc(sizeof(CSRbCSR));

    C->LLrow = malloc((n+1)*sizeof(int));
    //C->LLcol = malloc(100*sizeof(int));

    csr2bcsr(row,col,n,nnz,b,A,1);

    SpBMM2(A->LLcol, A->LLrow, A->Bcol, A->Brow, A->n, A->nnz, A->nnzb,
          A->LLrow, A->LLcol, A->Brow, A->Bcol, A->n, A->nnz, A->nnzb,
          &C->LLcol, C->LLrow, C->Bcol, C->Brow, C->n, &C->nnz, b, 2);

}


void test_mask(int *col, int *row, int N, int b){

    int csize = 1000;
    int *Crow_serial = malloc((N+1)*sizeof(int));
    int *Ccol_serial = malloc(csize*sizeof(int));

    BSpMM_masked(col,row,N,row,col,N,col,row,N,&Ccol_serial,Crow_serial,N,&csize,2);

}


void benchmark(int *col, int *row, int N, int b){

    struct timespec bS,bE,mBS,mBE,mS,mE;
    double bT, mBT, mT;

    //printCSR(row,col,n,n,b);

    CSRbCSR *A = malloc(sizeof(CSRbCSR));
    CSRbCSR *C = malloc(sizeof(CSRbCSR));

    C->LLrow = malloc((N+1)*sizeof(int));
    //C->LLcol = malloc(100*sizeof(int));


    clock_gettime(CLOCK_MONOTONIC, &bS);

    csr2bcsr(row,col,N,row[N],b,A,0);

    clock_gettime(CLOCK_MONOTONIC, &bE);
    bT = timeConv(diff(bS,bE));

    printf("build done\n");


    clock_gettime(CLOCK_MONOTONIC, &mBS);

    SpBMM2(A->LLcol, A->LLrow, A->Bcol, A->Brow, A->n, A->nnz, A->nnzb, //TODO:add print as option
          A->LLrow, A->LLcol, A->Brow, A->Bcol, A->n, A->nnz, A->nnzb,
          &C->LLcol, C->LLrow, C->Bcol, C->Brow, C->n, &C->nnz, b, 0);
    
    clock_gettime(CLOCK_MONOTONIC, &mBE);
    mBT = timeConv(diff(mBS,mBE));

    free(A->LLcol);
    free(A->LLrow);
    free(A->Bcol);
    free(A->Brow);
    free(A);

    printf("blocked done\n");


    clock_gettime(CLOCK_MONOTONIC, &mS);

    int csize = 1000;
    int *Crow_serial = malloc((N+1)*sizeof(int));
    int *Ccol_serial = malloc(csize*sizeof(int));

    //BSpMM3(col,row,N,row,col,N,&Ccol_serial,Crow_serial,N,&csize,0);

    clock_gettime(CLOCK_MONOTONIC, &mE);
    mT = timeConv(diff(mS,mE));

    printf("serial done");

    

    bool same = compare_CSR(C->LLcol,C->LLrow,Ccol_serial,Crow_serial,N,0);

    free(Ccol_serial);
    free(Crow_serial);
    free(C->LLrow);
    free(C->LLcol);
    free(C);
    

    printf("\nmatrices same: %d\n",same);

    printf("Times:\n Build=%lf, blocked=%lf, serial=%lf\n",bT,mBT,mT);

}

void test3v4(int *col, int *row, int N, int b){
    struct timespec bS,bE,mBS,mBE,mS,mE;
    double bT, mBT, mT;

    //printCSR(row,col,n,n,b);

    CSRbCSR *A = malloc(sizeof(CSRbCSR));
    CSRbCSR *C = malloc(sizeof(CSRbCSR));
    CSRbCSR *C2 = malloc(sizeof(CSRbCSR));
    

    C->LLrow = malloc((N+1)*sizeof(int));
    //C->LLcol = malloc(100*sizeof(int));
    C2->LLrow = malloc((N+1)*sizeof(int));


    clock_gettime(CLOCK_MONOTONIC, &bS);

    csr2bcsr(row,col,N,row[N],b,A,0);

    clock_gettime(CLOCK_MONOTONIC, &bE);
    bT = timeConv(diff(bS,bE));

    printf("build done\n");


    clock_gettime(CLOCK_MONOTONIC, &mBS);

    SpBMM2(A->LLcol, A->LLrow, A->Bcol, A->Brow, A->n, A->nnz, A->nnzb, //TODO:add print as option
          A->LLrow, A->LLcol, A->Brow, A->Bcol, A->n, A->nnz, A->nnzb,
          &C->LLcol, C->LLrow, C->Bcol, C->Brow, C->n, &C->nnz, b, 0);
    
    clock_gettime(CLOCK_MONOTONIC, &mBE);
    mBT = timeConv(diff(mBS,mBE));

    

    printf("blocked1 done\n");


    clock_gettime(CLOCK_MONOTONIC, &mS);

    SpBMM3(A->LLcol, A->LLrow, A->Bcol, A->Brow, A->n, A->nnz, A->nnzb, //TODO:add print as option
          A->LLcol, A->LLrow, A->Bcol, A->Brow, A->n, A->nnz, A->nnzb,
          &C2->LLcol, C2->LLrow, C2->Bcol, C2->Brow, C2->n, &C2->nnz, b, 0);


    clock_gettime(CLOCK_MONOTONIC, &mE);
    mT = timeConv(diff(mS,mE));

    printf("blocked2 done");

    

    bool same = compare_CSR(C->LLcol,C->LLrow,C2->LLcol,C2->LLrow,N,0);

    free(A->LLcol);
    free(A->LLrow);
    free(A->Bcol);
    free(A->Brow);
    free(A);

    free(C->LLrow);
    free(C->LLcol);
    free(C2->LLrow);
    free(C2->LLcol);
    free(C);
    

    printf("\nmatrices same: %d\n",same);

    printf("Times:\n Build=%lf, blocked1=%lf, blocked2=%lf\n",bT,mBT,mT);
}

void test10(int argc, char const *argv[]){
    
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nnz;
    uint32_t *I, *J;
    uint32_t *c3,*values;
    uint32_t isOneBased=0;

    int b = atoi(argv[2]);

    if ((f = fopen(argv[1], "r")) == NULL) 
        exit(1);
    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nnz)) !=0)
        exit(1);

    I = malloc(nnz *2 * sizeof(uint32_t));
    J = malloc(nnz *2 * sizeof(uint32_t));

    for (int i=0; i<nnz; i++)
    {
        fscanf(f, "%u %u\n", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
        //I[i+nnz]=J[i];
        //J[i+nnz]=I[i];
    }
    if (f !=stdin) fclose(f);

    uint32_t * row  = malloc((M + 1) * sizeof(uint32_t));   
    uint32_t * col  = malloc(nnz  * sizeof(uint32_t));

    coo2csc(col, row, I,J,nnz,M,isOneBased);

    free(I);
    free(J);

    //benchmark(col,row,M,b);
    //printCSR(row,col,M,M,M);
    //test_mask(col,row,M,b);
    test3v4(col,row,M,b);

    free(row);
    free(col);




}

void test9(){
    int col[] = {1,7,2,6,2,3,5,7,8,0,1,4,3};
    int row[] = {0,2,4,5,6,9,9,10,12,13};
    int n = 9;
    int nnz = 13;
    int b=3;

    benchmark(col,row,n,b);

}

void test8_4(){ 
    int col[] = {0,5,7,1,5,0,0,3,2,5,6};
    int row[] = {0,1,3,3,5,5,6,8,11};
    int n=8,b=4;

    benchmark(col,row,n,b);
}

void test_f2c(){
    bool C[] = {1,1,0,1,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0};
    int n=7, m=7, nnz=15;
    int *row =malloc((n+1)*sizeof(int));
    int *col = malloc(nnz*sizeof(int));

    int *row_t =malloc((n+1)*sizeof(int));
    int *col_t = malloc(nnz*sizeof(int));

    full2CSR(C,n,m,col,row);
    printCSR(row,col,n,m,n);

    full2CSR_T(C,n,m,col_t,row_t);
    printCSR(row_t,col_t,n,m,n);



    printf("\n");

    int init_size = 100;
    int *Crow = malloc((n+1)*sizeof(int));
    int *Ccol = malloc(init_size*sizeof(int));

    int *Crow_o = malloc((n+1)*sizeof(int));
    int *Ccol_o = malloc(init_size*sizeof(int));


    BSpMM3(col,row,n,row_t,col_t,n,&Ccol_o,Crow_o,n,&init_size,0);
    printCSR(Crow_o,Ccol_o,n,n,n);

    BSpMM4(col,row,n,col,row,n,&Ccol,Crow,n,&init_size,0);
    printCSR(Crow,Ccol,n,n,n);



}


void test_mb(){
    int th = 100000;
    int size = th*th;
    bool *hi = malloc(size*sizeof(bool));
    hi[0] = 0;
    for(int i=1; i<size; i++){
        hi[i] = hi[i-1] + 1;
    }
    sleep(5);
    free(hi);
}
int main(int argc, char const *argv[])
{
    //test10(argc,argv);
    //test8_4();
    //test_f2c();
    test_mb();
    return 0;
}
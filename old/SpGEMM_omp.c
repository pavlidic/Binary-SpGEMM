#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include "utils.h"
#include "ins.h"

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

bool SpGEMM_d(int  *Acol, int *Arow, int An, 
              int  *Bcol, int *Brow, int Bm,
              int **Ccol, int *Crow, int *Csize)
{
    //printCSR(Arow,Acol,An,An,An);
    int nnzcum=0;
    bool *xb = calloc(An,sizeof(bool));
    for(int i=0; i<An; i++){
        int nnzpv = nnzcum;//nnz of previous row;
        Crow[i] = nnzcum;
        if(nnzcum + An > *Csize){
            *Csize += MAX(An, *Csize/4);
            *Ccol = realloc(*Ccol,*Csize*sizeof(int));
        }

        for(int jj=Arow[i]; jj<Arow[i+1]; jj++){
            int j = Acol[jj];

            for(int kp=Brow[j]; kp<Brow[j+1]; kp++){
                int k = Bcol[kp];
                if(!xb[k]){
                    xb[k] = true;
                    (*Ccol)[nnzcum] = k;
                    nnzcum++;
                }

            }

        }
        if(nnzcum > nnzpv){
            quickSort(*Ccol,nnzpv,nnzcum-1);
            for(int p=nnzpv; p<nnzcum; p++){
                xb[ (*Ccol)[p] ] = false;
            }
        }

    }
    Crow[An] = nnzcum;

    free(xb);
    return Crow[An];
}


//SpGEMM for C(start_row:end_row,:)
bool SpGEMM_bigslice(int  *Acol, int *Arow, int An, 
                     int  *Bcol, int *Brow, int Bm,
                     int **Ccol, int *Crow, int *Csize,
                     int start_row, int end_row)
{
    //printCSR(Arow,Acol,An,An,An);
    int nnzcum=0;
    bool *xb = calloc(Bm,sizeof(bool));
    int ip=0;
    for(int i=start_row; i<end_row; i++){
        int nnzpv = nnzcum;//nnz of previous row;
        Crow[ip++] = nnzcum;

        if(nnzcum + Bm > *Csize){
            *Csize += MAX(Bm, *Csize/4);
            *Ccol = realloc(*Ccol,*Csize*sizeof(int));
        }

        for(int jj=Arow[i]; jj<Arow[i+1]; jj++){
            int j = Acol[jj];

            for(int kp=Brow[j]; kp<Brow[j+1]; kp++){
                int k = Bcol[kp];
                if(!xb[k]){
                    xb[k] = true;
                    (*Ccol)[nnzcum] = k;
                    nnzcum++;
                }
            }

        }
        if(nnzcum > nnzpv){
            quickSort(*Ccol,nnzpv,nnzcum-1);
            for(int p=nnzpv; p<nnzcum; p++){
                xb[ (*Ccol)[p] ] = false;
            }
        }

    }
    Crow[ip] = nnzcum;
    
    free(xb);
    /* #pragma omp critical
    {
    printCSR(Crow,*Ccol,end_row-start_row,Bm,Bm);
    } */
    return Crow[ip];
}


void SpGEMM_sliced(int  *Acol, int *Arow, int An, 
                   int  *Bcol, int *Brow, int Bm,
                   int **Ccol, int *Crow,
                   int slices)
{
    int threads = slices;//omp_get_num_threads();
    printf("threads=%d\n",threads);

    int tBlock = An/threads;

    int  *Ccol_sizes  = malloc(threads*sizeof(int ));
    int **Ccol_tBlock = malloc(threads*sizeof(int*));
    int **Crow_tBlock = malloc(threads*sizeof(int*));

    

    uint32_t i, init_size = Bm;

    for(i=0; i<threads; i++){
        Ccol_sizes[i]  = init_size;
        Ccol_tBlock[i] = malloc(init_size*sizeof(int));
        Crow_tBlock[i] = malloc((tBlock+1)*sizeof(int));
    }


    for( i=0; i<threads; i++ ){
        SpGEMM_bigslice(Acol,Arow,An,
                        Bcol,Brow,Bm,
                        &Ccol_tBlock[i],Crow_tBlock[i],&Ccol_sizes[i],
                        i*tBlock, (i+1)*tBlock);
        
    }


    int nnz =0;
    for(i=0; i<threads; i++){
        nnz += Crow_tBlock[i][tBlock];
    }

    *Ccol = malloc(nnz*sizeof(int));

    Crow[0] = 0;
    int nnzcum = 0;
    for(i=0; i<threads; i++){
        memcpy(&(*Ccol)[nnzcum], Ccol_tBlock[i], Crow_tBlock[i][tBlock]*sizeof(int));
        memcpy(&Crow[i*tBlock+1], &Crow_tBlock[i][1], tBlock*sizeof(int));
        nnzcum += Crow_tBlock[i][tBlock];
    }

    int row_sum = Crow[tBlock];
    for(i=1; i<threads; i++){
        for(int j=1; j<tBlock+1; j++){
            Crow[i*tBlock + j] += row_sum;
        }
        row_sum = Crow[(i+1)*tBlock];
    }

    

}

void SpGEMM_omp(int  *Acol, int *Arow, int An, 
                int  *Bcol, int *Brow, int Bm,
                int **Ccol, int *Crow,
                int slices)
{
    int threads = omp_get_num_threads();
    printf("threads=%d\n",threads);

    

    int tBlock = An/slices;

    int  *Ccol_sizes  = malloc(slices*sizeof(int ));
    int **Ccol_tBlock = malloc(slices*sizeof(int*));
    int **Crow_tBlock = malloc(slices*sizeof(int*));

    

    uint32_t i, init_size = Bm;

    for(i=0; i<slices; i++){
        Ccol_sizes[i]  = init_size;
        Ccol_tBlock[i] = malloc(init_size*sizeof(int));
        Crow_tBlock[i] = malloc((tBlock+1)*sizeof(int));
    }

    #pragma omp parallel private(i)
    {

        #pragma omp for schedule(static) nowait
        for( i=0; i<slices; i++ ){
            //printf("hi from %d\n",omp_get_thread_num());
            SpGEMM_bigslice(Acol,Arow,An,
                            Bcol,Brow,Bm,
                            &Ccol_tBlock[i],Crow_tBlock[i],&Ccol_sizes[i],
                            i*tBlock, (i+1)*tBlock);

        }
    }

    int nnz =0;
    for(i=0; i<slices; i++){
        nnz += Crow_tBlock[i][tBlock];
    }

    *Ccol = malloc(nnz*sizeof(int));

    Crow[0] = 0;
    int nnzcum = 0;
    for(i=0; i<slices; i++){
        memcpy(&(*Ccol)[nnzcum], Ccol_tBlock[i], Crow_tBlock[i][tBlock]*sizeof(int));
        memcpy(&Crow[i*tBlock+1], &Crow_tBlock[i][1], tBlock*sizeof(int));
        nnzcum += Crow_tBlock[i][tBlock];
    }

    int row_sum = Crow[tBlock];
    for(i=1; i<slices; i++){
        for(int j=1; j<tBlock+1; j++){
            Crow[i*tBlock + j] += row_sum;
        }
        row_sum = Crow[(i+1)*tBlock];
    }

    

}

void SpGEMM_omp2(int  *Acol, int *Arow, int An, 
                int  *Bcol, int *Brow, int Bm,
                int **Ccol, int *Crow,
                int tBlock)
{

    

    int slices = An/tBlock;
    int left = An%tBlock;

    int  *Ccol_sizes  = malloc(slices*sizeof(int ));
    int **Ccol_tBlock = malloc(slices*sizeof(int*));
    int **Crow_tBlock = malloc(slices*sizeof(int*));

    

    uint32_t i, init_size = 2*Bm;

    for(i=0; i<slices; i++){
        Ccol_sizes[i]  = init_size;
        Ccol_tBlock[i] = malloc(init_size*sizeof(int));
        Crow_tBlock[i] = malloc((tBlock+1)*sizeof(int));
    }

    tic;
    #pragma omp parallel private(i)
    {

        #pragma omp for schedule(dynamic) nowait
        for( i=0; i<slices; i++ ){
            //printf("hi from %d\n",omp_get_thread_num());
            SpGEMM_bigslice(Acol,Arow,An,
                            Bcol,Brow,Bm,
                            &Ccol_tBlock[i],Crow_tBlock[i],&Ccol_sizes[i],
                            i*tBlock, (i+1)*tBlock);

        }
    }
    printf("mult took: %lf\n",toc);
    tic;

    int nnz =0;
    for(i=0; i<slices; i++){
        nnz += Crow_tBlock[i][tBlock];
    }

    *Ccol = malloc(nnz*sizeof(int));

    Crow[0] = 0;
    int nnzcum = 0;
    for(i=0; i<slices; i++){
        memcpy(&(*Ccol)[nnzcum], Ccol_tBlock[i], Crow_tBlock[i][tBlock]*sizeof(int));
        memcpy(&Crow[i*tBlock+1], &Crow_tBlock[i][1], tBlock*sizeof(int));
        nnzcum += Crow_tBlock[i][tBlock];

        free(Ccol_tBlock[i]);
        free(Crow_tBlock[i]);
    }
    free(Ccol_sizes);

    int row_sum = Crow[tBlock];
    for(i=1; i<slices; i++){
        for(int j=1; j<tBlock+1; j++){
            Crow[i*tBlock + j] += row_sum;
        }
        row_sum = Crow[(i+1)*tBlock];
    }
    printf("rest took: %lf\n",toc);

    

}

void SpGEMM_omp3(int  *Acol, int *Arow, int An, 
                 int  *Bcol, int *Brow, int Bm,
                 int **Ccol, int *Crow,
                 int tBlock)
{

    

    int slices = An/tBlock;
    int left = An%tBlock;

    int  *Ccol_sizes  = malloc(slices*sizeof(int ));
    int **Ccol_tBlock = malloc(slices*sizeof(int*));
    int **Crow_tBlock = malloc(slices*sizeof(int*));

    

    uint32_t i, init_size = Bm;

    for(i=0; i<slices; i++){
        Ccol_sizes[i]  = init_size;
        Ccol_tBlock[i] = malloc(init_size*sizeof(int));
        Crow_tBlock[i] = malloc((tBlock+1)*sizeof(int));
    }

    int nnz=0;
    int *rowsum = malloc(slices*sizeof(int));
    int *slicesize = malloc(slices*sizeof(int));
    //tic;
    #pragma omp parallel default(shared) private(i)
    {
        int tid = omp_get_thread_num();

        #pragma omp for schedule(dynamic) reduction(+:nnz)
        for( i=0; i<slices; i++ ){
            SpGEMM_bigslice(Acol,Arow,An,
                            Bcol,Brow,Bm,
                            &Ccol_tBlock[i],Crow_tBlock[i],&Ccol_sizes[i],
                            i*tBlock, (i+1)*tBlock);
            rowsum[i] = Crow_tBlock[i][tBlock];
            nnz += Crow_tBlock[i][tBlock];
            slicesize[i] = Crow_tBlock[i][tBlock];
        }
        if(tid==0){
            for(i=1; i<slices; i++){
                rowsum[i] += rowsum[i-1];
                //printf("%d ",rowsum[i]);
            }
        }
        #pragma omp barrier
        if(tid!=0 && tid<slices){
            for(i=0; i<tBlock;i++){
                Crow_tBlock[tid][i+1] += rowsum[tid-1];
            }
        }   

    }


    *Ccol = malloc(nnz*sizeof(int));

    Crow[0] = 0;

    memcpy(&(*Ccol)[0], Ccol_tBlock[0], slicesize[0]*sizeof(int));
    memcpy(&Crow[1], &Crow_tBlock[0][1], tBlock*sizeof(int));

    free(Ccol_tBlock[0]);
    free(Crow_tBlock[0]);

    for(i=1; i<slices; i++){
        memcpy(&(*Ccol)[rowsum[i-1]], Ccol_tBlock[i], slicesize[i]*sizeof(int));
        memcpy(&Crow[i*tBlock+1], &Crow_tBlock[i][1], tBlock*sizeof(int));

        free(Ccol_tBlock[i]);
        free(Crow_tBlock[i]);
    }
    free(Ccol_sizes);
    free(Ccol_tBlock);
    free(Crow_tBlock);

}



bool SpGEMM_d_masked( int  *Acol, int *Arow, int An, 
                        int  *Bcol, int *Brow, int Bm,
                        int  *Fcol, int *Frow,
                        int **Ccol, int *Crow, int *Csize)//output

{
    //printCSR(Arow,Acol,An,An,An);
    int nnzcum=0;
    bool *xb = calloc(An,sizeof(bool));

    for(int i=0; i<An; i++){
        int nnzpv = nnzcum;//nnz of previous row;
        Crow[i] = nnzcum;
        if(nnzcum + An > *Csize){
            *Csize += MAX(An, *Csize/4);
            *Ccol = realloc(*Ccol,*Csize*sizeof(int));
        }

        //add row mask
        for(int jj=Frow[i]; jj<Frow[i+1]; jj++){
            xb[ Fcol[jj] ] = true;
        }


        //add new row items
        for(int jj=Arow[i]; jj<Arow[i+1]; jj++){
            int j = Acol[jj];

            for(int kp=Brow[j]; kp<Brow[j+1]; kp++){
                int k = Bcol[kp];
                if(!xb[k]){
                    xb[k] = true;
                    (*Ccol)[nnzcum] = k;
                    nnzcum++;
                }
            }

        }

        if(nnzcum > nnzpv){
            quickSort(*Ccol,nnzpv,nnzcum-1);
            for(int p=nnzpv; p<nnzcum; p++){
                xb[ (*Ccol)[p] ] = false;
            }
        }

        //turn of row mask bits;
        for(int jj=Frow[i]; jj<Frow[i+1]; jj++){
            xb[ Fcol[jj] ] = false;
        }

    }
    Crow[An] = nnzcum;

    free(xb);
    return nnzcum;
}




int test_slice(int argc, char const *argv[]){
    struct timespec bS,bE,nbS,nbE;
    int s = atoi(argv[2]);
    int e = atoi(argv[3]);
    int b = 3;
    
    uint32_t *Arow,*Acol;
    int An,Am,Annz;
    readCOO(argv[1],&Arow,&Acol,&An,&Am,&Annz);
    printf("n=%d, nnz=%d\n",An,Annz);
    //printCSR(Arow,Acol,An,Am,b);




    int Csize = 1000;
    int *nCrow = calloc((An+1),sizeof(int));
    int *nCcol = malloc(Csize*sizeof(int));
    clock_gettime(CLOCK_MONOTONIC, &nbS);
    SpGEMM_d(Acol,Arow,An,Acol,Arow,An,&nCcol,nCrow,&Csize);
    clock_gettime(CLOCK_MONOTONIC, &nbE);
    printf("whole: %lf\n",timeConv(diff(nbS,nbE)));
    printf("nnz=%d\n",nCrow[An]);
    printCSR(nCrow,nCcol,An,An,b);

    free(nCcol);
    free(nCrow);

    Csize = 1000;
    nCrow = calloc((An+1),sizeof(int));
    nCcol = malloc(Csize*sizeof(int));
    clock_gettime(CLOCK_MONOTONIC, &nbS);
    SpGEMM_bigslice(Acol,Arow,An,Acol,Arow,An,&nCcol,nCrow,&Csize,s,e);
    clock_gettime(CLOCK_MONOTONIC, &nbE);
    printf("slice: %lf\n",timeConv(diff(nbS,nbE)));
    printf("nnz=%d\n",nCrow[An]);
    printCSR(nCrow,nCcol,e-s,An,b);


}

int test_omp(int argc, char const *argv[]){
    struct timespec bS,bE,nbS,nbE;
    int slices = atoi(argv[2]);
    printf("nt=%d\n",slices);
    int b = 3;

    omp_set_num_threads(slices);
    
    uint32_t *Arow,*Acol;
    int An,Am,Annz;
    readCOO(argv[1],&Arow,&Acol,&An,&Am,&Annz);
    printf("n=%d, nnz=%d\n",An,Annz);
    //printCSR(Arow,Acol,An,Am,b);




    int Csize = 1000;
    int *nCrow = calloc((An+1),sizeof(int));
    int *nCcol = malloc(Csize*sizeof(int));
    clock_gettime(CLOCK_MONOTONIC, &nbS);
    SpGEMM_omp(Acol,Arow,An,Acol,Arow,An,&nCcol,nCrow,slices);
    clock_gettime(CLOCK_MONOTONIC, &nbE);
    printf("whole: %lf\n",timeConv(diff(nbS,nbE)));
    printf("nnz=%d\n",nCrow[An]);
    //printCSR(nCrow,nCcol,An,An,b);

    free(nCcol);
    free(nCrow);

}


int test_omp2(int argc, char const *argv[]){
    struct timespec bS,bE,nbS,nbE;
    int tBlock = atoi(argv[2]);
    int threads = atoi(argv[3]);
    printf("threads=%d, tBlock=%d\n",threads,tBlock);
    int b = 3;

    omp_set_num_threads(threads);
    
    uint32_t *Arow,*Acol;
    int An,Am,Annz;
    readCOO(argv[1],&Arow,&Acol,&An,&Am,&Annz);
    printf("n=%d, nnz=%d\n",An,Annz);
    //printCSR(Arow,Acol,An,Am,b);




    int Csize = 1000;
    int *nCrow = calloc((An+1),sizeof(int));
    int *nCcol = malloc(Csize*sizeof(int));
    clock_gettime(CLOCK_MONOTONIC, &nbS);
    SpGEMM_omp2(Acol,Arow,An,Acol,Arow,An,&nCcol,nCrow,tBlock);
    clock_gettime(CLOCK_MONOTONIC, &nbE);
    printf("whole: %lf\n",timeConv(diff(nbS,nbE)));
    printf("nnz=%d\n",nCrow[An]);
    //printCSR(nCrow,nCcol,An,An,b);

    free(nCcol);
    free(nCrow);

}

int test_sliced(int argc, char const *argv[]){
    struct timespec bS,bE,nbS,nbE;
    int slices = atoi(argv[2]);

    int b = 3;

    
    uint32_t *Arow,*Acol;
    int An,Am,Annz;
    readCOO(argv[1],&Arow,&Acol,&An,&Am,&Annz);
    printf("n=%d, nnz=%d\n",An,Annz);
    //printCSR(Arow,Acol,An,Am,b);




    int Csize = 1000;
    int *nCrow = calloc((An+1),sizeof(int));
    int *nCcol = malloc(Csize*sizeof(int));
    clock_gettime(CLOCK_MONOTONIC, &nbS);
    SpGEMM_sliced(Acol,Arow,An,Acol,Arow,An,&nCcol,nCrow,slices);
    clock_gettime(CLOCK_MONOTONIC, &nbE);
    printf("whole: %lf\n",timeConv(diff(nbS,nbE)));
    printf("nnz=%d\n",nCrow[An]);
    printCSR(nCrow,nCcol,An,An,b);

    free(nCcol);
    free(nCrow);

}

int test_randcsr(int argc, char const *argv[]){
    srand(time(NULL));
    int n = atoi(argv[1]);
    int d = atoi(argv[2]);


    int *Ccol,*Crow;

    tic;
    rndCSR(&Crow,&Ccol,n,n,d);
    //printCSR(Crow,Ccol,n,n,10);
    printf("nnz=%d,time=%lf\n",Crow[n],toc);

}
int main(int argc, char const *argv[])
{
    //printf("%d\n",argc);
    if(argc<2){
        printf("usage: ./BSpGEMM matrix block-size\n");
        exit(1);
    }
    //test_slice(argc,argv);
    test_omp2(argc,argv);
    //test_randcsr(argc,argv);
    //test_sliced(argc,argv);
    return 0;
}


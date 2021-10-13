#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include "utils.h"
#include "ins.h"

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define isroot if(rank==0)

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

    

    uint32_t i, init_size = Bm;

    for(i=0; i<slices; i++){
        Ccol_sizes[i]  = init_size;
        Ccol_tBlock[i] = malloc(init_size*sizeof(int));
        Crow_tBlock[i] = malloc((tBlock+1)*sizeof(int));
    }

    //tic;
    #pragma omp parallel private(i)
    {
        //int tid=omp_get_thread_num();tic;

        #pragma omp for schedule(static) nowait
        for( i=0; i<slices; i++ ){
            //printf("hi from %d\n",omp_get_thread_num());
            SpGEMM_bigslice(Acol,Arow,An,
                            Bcol,Brow,Bm,
                            &Ccol_tBlock[i],Crow_tBlock[i],&Ccol_sizes[i],
                            i*tBlock, (i+1)*tBlock);

        }
        //printf("thread:%d, time:%lf\n",tid,toc);
    }
    //printf("mult took: %lf\n",toc);
    //tic;

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
    free(Ccol_tBlock);
    free(Crow_tBlock);

    int row_sum = Crow[tBlock];
    for(i=1; i<slices; i++){
        for(int j=1; j<tBlock+1; j++){
            Crow[i*tBlock + j] += row_sum;
        }
        row_sum = Crow[(i+1)*tBlock];
    }
    //printf("rest took: %lf\n",toc);
    //printCSR(Crow,*Ccol,An,Bm,Bm);

    

}

void SpGEMM_mpi(int  *Acol, int *Arow, int An, 
                int  *Bcol, int *Brow, int Bm,
                int **Ccol, int *Crow,
                int tBlock)
{
    int numtasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int tasksize = An/numtasks;     //number of rows an mpi task will have to compute

    int *Ccol_slice;
    int *Crow_slice = malloc((tasksize+1)*sizeof(int));

    //tic;
    SpGEMM_omp2(Acol,&Arow[rank*tasksize],tasksize,
                Bcol,Brow,Bm,
                &Ccol_slice,Crow_slice,
                tBlock);//TODO: every task on their own matrices, then gather like below
    //printf("toc:%lf\n",toc);
    //printf("hifrom: %d\n",rank);
    //printCSR(Crow_slice,Ccol_slice,tasksize,Bm,Bm);

    
    //find final nnz
    int all_nnz=0;
    MPI_Reduce(&Crow_slice[tasksize],&all_nnz,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    
    
    //find nnz per task
    int *reccounts;
    isroot{
        reccounts = malloc(numtasks*sizeof(int));
    }
    MPI_Gather(&Crow_slice[tasksize],1,MPI_INT,reccounts,1,MPI_INT,0,MPI_COMM_WORLD);

    //compute displacement array
    int *disps;
    isroot{
        disps = malloc(numtasks*sizeof(int));
        disps[0] = 0;
        for(int t=1; t<numtasks; t++){
            disps[t] = disps[t-1] + reccounts[t-1];
        }
    }



    //gather results of each task
    isroot{
        //printf("allnnz=%d\n",all_nnz);
        *Ccol = malloc(all_nnz*sizeof(int));
        //Crow  = malloc((An+1)*sizeof(int));
        Crow[0] = 0;
    }
    MPI_Gatherv(Ccol_slice,Crow_slice[tasksize],MPI_INT,*Ccol,reccounts,disps,MPI_INT,0,MPI_COMM_WORLD);//gather Ccol
    MPI_Gather(&Crow_slice[1],tasksize,MPI_INT,&Crow[1],tasksize,MPI_INT,0,MPI_COMM_WORLD);             //gather Crow
    
    free(Crow_slice);
    free(Ccol_slice);


    //fix Crow from representing each slice to representing the whole matrix
    isroot{
        free(disps);
        free(reccounts);

        int row_sum = Crow[tasksize];
        for(int i=1; i<numtasks; i++){
            for(int j=1; j<tasksize+1; j++){
                Crow[i*tasksize + j] += row_sum;
            }
            row_sum = Crow[(i+1)*tasksize];
        }
        //printCSR(Crow,*Ccol,An,Bm,Bm);
    }

    /* isroot{
        for(int i=0; i<numtasks; i++){
            for(int j=1; j<tasksize+1; j++){
                printf("%d ",Crow[i*tasksize + j]);
            }
            printf("\n");
        }
    } */
    

    

}

void SpGEMM_mpi_omp(int  *Acol, int *Arow, int An, 
                    int  *Bcol, int *Brow, int Bm,
                    int **Ccol, int *Crow,
                    int tBlock)
{
    int numtasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int tasksize = An/numtasks; //number of rows an mpi task will have to compute

    int slices = tasksize/tBlock; //number of blocks of rows every task will have to compute


    int  *Ccol_sizes  = malloc(slices*sizeof(int ));
    int **Ccol_tBlock = malloc(slices*sizeof(int*));
    int **Crow_tBlock = malloc(slices*sizeof(int*));
    

    uint32_t i, init_size = Bm;

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
            SpGEMM_bigslice(Acol,&Arow[rank*tasksize],slices,
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
    }

    int row_sum = Crow[tBlock];
    for(i=1; i<slices; i++){
        for(int j=1; j<tBlock+1; j++){
            Crow[i*tBlock + j] += row_sum;
        }
        row_sum = Crow[(i+1)*tBlock];
    }
    printf("rest took: %lf\n",toc);
    //printCSR(Crow,*Ccol,tasksize,Bm,Bm);


    int all_nnz=0;
    MPI_Reduce(&Crow[tasksize],&all_nnz,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(rank==0){
        printf("allnnz=%d\n",all_nnz);
    }

    int *tempCcol = malloc(all_nnz*sizeof(int));
    int *tempCrow = malloc((An+1)*sizeof(int));
    tempCrow[0]=0;

    MPI_Gather(*Ccol,tasksize,MPI_INT,tempCcol,tasksize,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&Crow[1],tasksize,MPI_INT,&tempCrow[1],tasksize,MPI_INT,0,MPI_COMM_WORLD);


    if(rank==0){
        printf("tasksize=%d\n",tasksize);
        

        int temp_row_sum = tempCrow[tasksize];
        for(i=1; i<slices; i++){
            for(int j=1; j<tasksize+1; j++){
                //printf("");
                tempCrow[i*tasksize + j] += temp_row_sum;
            }
            temp_row_sum = tempCrow[(i+1)*tasksize];
        }
        printCSR(tempCrow,tempCcol,An,Bm,Bm);
    }
    

    

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

int test_mpi(int argc, char const *argv[]){
    
    int numtasks, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    int tBlock  = atoi(argv[2]);
    int threads = atoi(argv[3]);
    int times   = atoi(argv[4]);

    omp_set_num_threads(threads);
    
    uint32_t *Arow,*Acol;
    int An,Am,Annz;
    readCOO(argv[1],&Arow,&Acol,&An,&Am,&Annz);
    //printCSR(Arow,Acol,An,Am,b);
    

    int *nCrow = calloc((An+1),sizeof(int));
    int *nCcol;

    
    double *alltimes = malloc(times*sizeof(double));
    double timesum = 0;

    for(int i=0; i<times; i++){
        MPI_Barrier(MPI_COMM_WORLD);
        tic;

        SpGEMM_mpi(Acol,Arow,An,Acol,Arow,An,&nCcol,nCrow,tBlock);

        alltimes[i] = toc;
        timesum += alltimes[i];

        isroot free(nCcol);
    }

    double mean = timesum/times;
    quickSortD(alltimes,0,times-1);
    double median  = alltimes[(times-1)/2];
    double fastest = alltimes[0];
    /* isroot{
        for(int i=0; i<times; i++){
            printf("i %d, all %lf \n",i,alltimes[i]);
        }
    } */    
    //printf("%d: %lf\n",rank,total_time);
    isroot{
        //printCSR(nCrow,nCcol,An,An,b);
        printf("%d,%d,%d,%d,%d,%d,%lf,%lf,%lf\n",numtasks,threads,tBlock,An,Annz,nCrow[An],mean,median,fastest);
    }
    free(alltimes);
    free(Acol);
    free(Arow);

    free(nCrow);

}

/* int test_mpi2(int argc, char const *argv[]){
    
    int numtasks, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double totaltime;


    int tBlock  = atoi(argv[2]);
    int threads = atoi(argv[3]);

    omp_set_num_threads(threads);
    
    uint32_t *Arow,*Acol;
    int An,Am,Annz;
    readCOO(argv[1],&Arow,&Acol,&An,&Am,&Annz);
    //printCSR(Arow,Acol,An,Am,b);
    

    int *nCrow = calloc((An+1),sizeof(int));
    int *nCcol;

    



    MPI_Barrier(MPI_COMM_WORLD);
    tic;

    SpGEMM_mpi(Acol,Arow,An,Acol,Arow,An,&nCcol,nCrow,tBlock);

    totaltime = toc;


    //printf("%d: %lf\n",rank,total_time);
    isroot{
        //printCSR(nCrow,nCcol,An,An,b);
        printf("%d,%d,%d,%d,%d,%d,%lf,\n",numtasks,threads,tBlock,An,Annz,nCrow[An],totaltime);
    }

    
    free(nCrow);

} */


int main(int argc, char const *argv[])
{
    int numtasks, rank, provided;

	//MPI_Init(&argc,&argv);
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    //printf("req:%d, prov:%d\n",MPI_THREAD_FUNNELED,provided);
    MPI_Query_thread(&provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //printf("%d\n",argc);
    if(argc<4){
        printf("usage: ./SpGEMM  path-to-matrix  threadslice_size  number_of_threads times_to_run\n");
        exit(1);
    }
    //test_slice(argc,argv);
    test_mpi(argc,argv);
    //test_sliced(argc,argv);

    MPI_Finalize();
    return 0;
}


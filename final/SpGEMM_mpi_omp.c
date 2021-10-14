#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include "utils.h"

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define isroot if(rank==0)


//SpGEMM for C(start_row:end_row,:)
void SpGEMM_bigslice(int  *Acol, int *Arow, int An, 
                     int  *Bcol, int *Brow, int Bm,
                     int **Ccol, int *Crow, int *Csize,
                     int start_row, int end_row)
{
    int nnzcum=0;                       // the sum of nnz at any given time
    bool *xb = calloc(Bm,sizeof(bool)); // a binary array of flags to not add elements more than one time
    int ip=0;                           // row pointer in Crow

    for(int i=start_row; i<end_row; i++){
        int nnzpv = nnzcum;             // nnz of previous row;
        Crow[ip++] = nnzcum;            // update Crow at the start of every row caclulation

        if(nnzcum + Bm > *Csize){       // make more space if there isnt at least enough for a whole row (the max elements that can be added)
            *Csize += MAX(Bm, *Csize/4);
            *Ccol = realloc(*Ccol,*Csize*sizeof(int));
        }

        for(int jj=Arow[i]; jj<Arow[i+1]; jj++){    // gustavson algorithm
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
        if(nnzcum > nnzpv){                         // if there were any added
        quickSort(*Ccol,nnzpv,nnzcum-1);            // sort the row since they could be out of order
            for(int p=nnzpv; p<nnzcum; p++){        // reset the flag bit array
                xb[ (*Ccol)[p] ] = false;
            }
        }

    }
    Crow[ip] = nnzcum;
    
    free(xb);

}


//OpenMP implementation of SpGEMM
//A and B are CSR matrices, C will also be CSR
//
//tBlock is the number of rows each thread will do
//
//needs n to be divisible by tBlock.
//
//making it work with non-divisible matrix sizes is easy to do
//but not required for the extent of this assignment
//since it will always be for test matrices with lengths multiples of 10
void SpGEMM_omp(int  *Acol, int *Arow, int An, 
                 int  *Bcol, int *Brow, int Bm,
                 int **Ccol, int *Crow,
                 int tBlock)
{

    int slices = An/tBlock; // number of slices it will be devided into

    // result will be placed into different CSR matrices and then joined together at the end
    int  *Ccol_sizes  = malloc(slices*sizeof(int ));    // saves the Ccol size of each CSR matrix
    int **Ccol_tBlock = malloc(slices*sizeof(int*));    // Ccol for every slice
    int **Crow_tBlock = malloc(slices*sizeof(int*));    // Crow for every slice


    uint32_t i, init_size = Bm; // initial size of every Ccol is Bm

    // initialize the matrices
    for(i=0; i<slices; i++){
        Ccol_sizes[i]  = init_size;
        Ccol_tBlock[i] = malloc(init_size*sizeof(int));
        Crow_tBlock[i] = malloc((tBlock+1)*sizeof(int));
    }

    // do every slice in parallel
    #pragma omp parallel private(i)
    {

        #pragma omp for schedule(static) nowait
        for( i=0; i<slices; i++ ){

            SpGEMM_bigslice(Acol,Arow,An,
                            Bcol,Brow,Bm,
                            &Ccol_tBlock[i],Crow_tBlock[i],&Ccol_sizes[i],
                            i*tBlock, (i+1)*tBlock);

        }

    }

    // calculate nnz for the final Ccol and allocate enough space at once
    int nnz =0;
    for(i=0; i<slices; i++){
        nnz += Crow_tBlock[i][tBlock];
    }
    *Ccol = malloc(nnz*sizeof(int));

    // construct the final result matrix by copying the slices
    Crow[0] = 0;
    int nnzcum = 0; //Ccol pointer for the next slice to copy to
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


    // fix Crow from representing each slice to representing the whole matrix
    int row_sum = Crow[tBlock];
    for(i=1; i<slices; i++){
        for(int j=1; j<tBlock+1; j++){
            Crow[i*tBlock + j] += row_sum;
        }
        row_sum = Crow[(i+1)*tBlock];
    }    

}

//MPI implementation of SpGEMM
//A and B are CSR matrices, C will also be CSR
//
//tBlock is the number of rows each thread will do
//
//needs n to be divisible by numtasks 
//
//making it work with non-divisible matrix sizes is easy to do
//but not required for the extent of this assignment
//since it will always be for test matrices with lengths multiples of 10
void SpGEMM_mpi(int  *Acol, int *Arow, int An, 
                int  *Bcol, int *Brow, int Bm,
                int **Ccol, int *Crow,
                int tBlock)
{
    // initialize mpi task variables
    int numtasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int tasksize = An/numtasks; // number of rows an mpi task will have to compute

    int *Ccol_slice;    // Ccol is alocated inside SpGEMM_omp exactly
    int *Crow_slice = malloc((tasksize+1)*sizeof(int));

    // call SpGEMM_omp only on the part each task has to do
    SpGEMM_omp(Acol,&Arow[rank*tasksize],tasksize,
                Bcol,Brow,Bm,
                &Ccol_slice,Crow_slice,
                tBlock);
    
    // find final nnz
    int all_nnz=0;
    MPI_Reduce(&Crow_slice[tasksize],&all_nnz,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    
    
    // find nnz per task needed for the Ccol gather
    int *reccounts;
    isroot{
        reccounts = malloc(numtasks*sizeof(int));
    }
    MPI_Gather(&Crow_slice[tasksize],1,MPI_INT,reccounts,1,MPI_INT,0,MPI_COMM_WORLD);

    // compute displacement array needed for the Ccol gather on the root task
    int *disps;
    isroot{
        disps = malloc(numtasks*sizeof(int));
        disps[0] = 0;
        for(int t=1; t<numtasks; t++){
            disps[t] = disps[t-1] + reccounts[t-1];
        }
    }

    // gather results of each task needed for the Ccol gather on the root task
    isroot{
        *Ccol = malloc(all_nnz*sizeof(int));
        Crow[0] = 0;
    }
    MPI_Gatherv(Ccol_slice,Crow_slice[tasksize],MPI_INT,*Ccol,reccounts,disps,MPI_INT,0,MPI_COMM_WORLD);//gather Ccol
    MPI_Gather(&Crow_slice[1],tasksize,MPI_INT,&Crow[1],tasksize,MPI_INT,0,MPI_COMM_WORLD);             //gather Crow
    
    free(Crow_slice);
    free(Ccol_slice);


    // fix Crow from representing each task to representing the whole matrix
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
        
    }

}


// serial implementation of masked SpGEMM
// can easily be parallelized exactly as above
// does C = F.*(A*B)
// all in CSR form
void SpGEMM_masked(int  *Acol, int *Arow, int An, 
                   int  *Bcol, int *Brow, int Bm,
                   int  *Fcol, int *Frow,
                   int **Ccol, int *Crow, int *Csize)

{

    int nnzcum=0;
    bool *xb = malloc(An*sizeof(bool));
    for(int i=0; i<An; i++) xb[i] = true;           // start at 1, then make false only the mask bits

    for(int i=0; i<An; i++){
        int nnzpv = nnzcum; //nnz of previous row;
        Crow[i] = nnzcum;
        if(nnzcum + An > *Csize){
            *Csize += MAX(An, *Csize/4);
            *Ccol = realloc(*Ccol,*Csize*sizeof(int));
        }

        //add row mask
        for(int jj=Frow[i]; jj<Frow[i+1]; jj++){    // this and the part below using F is the only difference
            xb[ Fcol[jj] ] = false;                 // to the normal algorithm.
        }                                           // It just turns the flag bits than can be added to 0


        // gustavson
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

        
        for(int jj=Frow[i]; jj<Frow[i+1]; jj++){   
            xb[ Fcol[jj] ] = true;             // reset the flag bits to 1
        }

    }
    Crow[An] = nnzcum;

    free(xb);

}




// test function
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

    int *nCrow = calloc((An+1),sizeof(int));
    int *nCcol;

    
    double *alltimes = malloc(times*sizeof(double));
    double timesum = 0;

    for(int i=0; i<times; i++){
        MPI_Barrier(MPI_COMM_WORLD);    // start the timer at the same time for all tasks
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
    isroot{
        // prints tasks,threads,total cpus,tbloc,matrix used,An,Annz,Cnnz,mean time,nedian time,fastest time
        printf("%d,%d,%d,%d,%s,%d,%d,%d,%lf,%lf,%lf\n",numtasks,threads,numtasks*threads,tBlock,argv[1],An,Annz,nCrow[An],mean,median,fastest);
    }
    free(alltimes);
    free(Acol);
    free(Arow);

    free(nCrow);

}

int main(int argc, char const *argv[])
{
    int numtasks, rank, provided;
    

    //initialize MPI
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    MPI_Query_thread(&provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc!=5){
        printf("usage: mpirun  -n  numtasks  SpGEMM_mpi_omp  path-to-matrix  threadslice_size  number_of_threads  times_to_run\n");
        exit(1);
    }
    test_mpi(argc,argv);
    

    MPI_Finalize();
    return 0;
}


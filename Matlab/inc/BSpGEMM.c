#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "BSpGEMM.h"
#include "utils.h"

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

void SpGEMM_mat(int *Acol, int *Arow, int An, 
                int *Bcol, int *Brow, int Bm,
                int *Ccol, int *Crow)
{
    //printCSR(Arow,Acol,An,An,An);


    int nnzcum=0;
    bool *xb = calloc(An,sizeof(bool));
    for(int i=0; i<An; i++){
        int nnzpv = nnzcum;//nnz of previous row;
        Crow[i] = nnzcum;

        for(int jj=Arow[i]; jj<Arow[i+1]; jj++){
            int j = Acol[jj];

            for(int kp=Brow[j]; kp<Brow[j+1]; kp++){
                int k = Bcol[kp];
                if(!xb[k]){
                    xb[k] = true;
                    Ccol[nnzcum] = k;
                    nnzcum++;
                }

            }

        }
        if(nnzcum > nnzpv){
            quickSort(Ccol,nnzpv,nnzcum-1);
            for(int p=nnzpv; p<nnzcum; p++){
                xb[ Ccol[p] ] = false;
            }
        }

    }
    Crow[An] = nnzcum;

    free(xb);
}


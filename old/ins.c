#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

//return 1 if found, 0 if not
bool binSearch(int array[], int first, int last, int search, int *pos){
    int middle = (first+last)/2;
    while(first <= last) {
        if (array[middle] < search)
            first = middle + 1;
        else if (array[middle] == search) {
            return 1;
        }
        else
            last = middle - 1;

        middle = (first + last)/2;
    }
    if(first > last){
        *pos = first;
        return 0;
    }
        
}

void insertCSR(int **Ccol, int *Crow, int item, int *Csize, int row, int n){
    int pos;
    bool found = binSearch(*Ccol, Crow[row], Crow[row+1] - 1, item, &pos);
    if(!found){
        int j = Crow[n];
        if(j==*Csize){ //if nnz = Csize (full)
            *Csize *= 2;
            *Ccol = realloc(*Ccol,*Csize*sizeof(int));
        }
        
        while(j>pos){
            (*Ccol)[j] = (*Ccol)[j-1];
            j--;
        }
        (*Ccol)[pos] = item;

        for(int i = row; i<n; i++){
            Crow[i+1]++;
        }

    }
}

void insertCSR_mv(int **Ccol, int *Crow, int item, int *Csize, int row, int n){
    int pos;
    const bool found = binSearch(*Ccol, Crow[row], Crow[row+1] - 1, item, &pos);
    if(!found){
        int j = Crow[n];
        if(j==*Csize){ //if nnz = Csize (full)
            *Csize *= 2;
            *Ccol = realloc(*Ccol,*Csize*sizeof(int));
        }
        
        memmove(&(*Ccol)[pos+1],&(*Ccol)[pos],(Crow[n]-pos)*sizeof(int));
        (*Ccol)[pos] = item;
        
        for(int i = row; i<n; i++){
            Crow[i+1]++;
        }

    }
}

bool insertCSR_mv_noadd(int **Ccol, int *Crow, int item, int *Csize, int row, int n, int added){
    int pos;
    const bool found = binSearch(*Ccol, Crow[row], Crow[row+1] + added - 1, item, &pos);
    if(!found){
        int j = Crow[n];
        if(j==*Csize){ //if nnz = Csize (full)
            *Csize *= 2;
            *Ccol = realloc(*Ccol,*Csize*sizeof(int));
        }
        
        memmove(&(*Ccol)[pos+1],&(*Ccol)[pos],(Crow[n]+added-pos)*sizeof(int));
        (*Ccol)[pos] = item;
        
        return 1;
    }
    return 0;
}


/* void test8_4(int argc, char const *argv[]){ 
    int col[] = {0,5,7,1,5,0,0,3,2,5,6};
    int row[] = {0,1,3,3,5,5,6,8,11};
    int n=8,nnz = 11;

    int Csize = nnz;
    int *Ccol = malloc(Csize*sizeof(int));
    int *Crow = malloc((n+1)*sizeof(int));

    for(int i=0; i<nnz; i++) Ccol[i]=col[i];
    for(int i=0; i<=n;  i++) Crow[i]=row[i];

    //printCSR(row,col,n,n,n);

    int item = atoi(argv[1]);
    int r = atoi(argv[2]);
    printf("\nitem to insert = %d in row =%d\n",item,r);
    insertCSR_mv(&Ccol,Crow,item,&Csize,r,n);
    //printCSR(Crow,Ccol,n,n,n);
}

void test_mv(){
    int col[] = {0,5,7,1,5,0,0,3,2,5,6};
    int n=8,nnz = 11;
    int Csize = nnz;
    int *Ccol = malloc(Csize*sizeof(int));
    for(int i=0; i<nnz; i++) Ccol[i]=col[i];
    for(int i=0; i<nnz; i++) printf("%d ",Ccol[i]);
    printf("\n");


    int pos = 3;

    memmove(&Ccol[pos+1],&Ccol[pos],(nnz-pos)*sizeof(int));
    Ccol[pos] = 33;
    for(int i=0; i<nnz+1; i++) printf("%d ",Ccol[i]);
    printf("\n");

}

void big_test(){
    int n=5;
    int Csize = 5;
    int *Ccol = malloc(Csize*sizeof(int));
    int *Crow = malloc((n+1)*sizeof(int));
    for(int i=0; i<=n;  i++) Crow[i]=0;
    printCSR(Crow,Ccol,n,n,n);

    int item, row;
    while(1){
        printf("\nnext to insert: ");
        scanf("%d %d",&item, &row);

        insertCSR(&Ccol,Crow,item,&Csize,row,n);
        printCSR(Crow,Ccol,n,n,n);
        printf("size =%d\n",Csize);
    }

}
 */

/* int main(int argc, char const *argv[])
{
    //test_mv();

    //test8_4(argc,argv);
    big_test();
    return 0;
} */


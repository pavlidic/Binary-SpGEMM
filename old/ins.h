#ifndef INS_H
#define INS_H

//binary search, returns 1 if found, 0 if not and puts position in pos
bool binSearch(const int array[], int first, int last, const int search, int *pos);

//inserts item in CSR matrix at (item,row) using loop
void insertCSR(int **Ccol, int *Crow, const int item, int *Csize, const int row, const int n);

//inserts item in CSR matrix at (item,row) using memmove
void insertCSR_mv(int **Ccol, int *Crow, const int item, int *Csize, const int row, const int n);

bool insertCSR_mv_noadd(int **Ccol, int *Crow, const int item, int *Csize, const int row, const int n, int added);

#endif
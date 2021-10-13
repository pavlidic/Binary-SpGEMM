function test_SpGEMM(n,d)
% n = collumn/row size
% d = nz per collumn/row

%generate random matrices
A=sprand(n,n,d/n);
B=sprand(n,n,d/n);


[An,Am] = size(A);
[Bn,Bm] = size(B);

%test if multiplication is possible
if(Am~=Bn)
    fprintf("\nSize missmatch!\n\n")
    return;
end

%do the matlab multiplication
tic;Cmat = A*B>0;mattime=toc;
Csize = nnz(Cmat);



%matlab-sparse to CSR
base = 0;
[~, Arow_ptr, Acol_ind] = sparse2csr(A, base);
[~, Brow_ptr, Bcol_ind] = sparse2csr(B, base);

%turn to int32
Arow_ptr = int32(Arow_ptr);
Brow_ptr = int32(Brow_ptr);
Acol_ind = int32(Acol_ind);
Bcol_ind = int32(Bcol_ind);


%generate C code if it doesnt already exist
codegen SpGEMM -args {Arow_ptr, Acol_ind, An, Brow_ptr, Bcol_ind, Bm, Csize} ./inc/BSpGEMM.c ./inc/utils.c


tic;
[Ccol_ind, Crow_ptr] = SpGEMM_mex(Arow_ptr, Acol_ind, An, Brow_ptr, Bcol_ind, Bm, Csize);
ctime=toc;

%generate array of 1s to use when rebuilding the sparse array C
val = ones(1,Csize);

%rebuild C (in C-language)
Cc = csr2sparse(val, double(Crow_ptr), double(Ccol_ind), An);


%test if equal
if(isequal(Cmat,Cc))
    fprintf("Matrices are the same!\n")
    fprintf("Matlab time: %d sec\n",mattime)
    fprintf("C      time: %d sec\n",ctime)
else
    fprintf("Matrices are different!\n")
end

end
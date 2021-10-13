function [Ccol_ind, Crow_ptr] = SpGEMM(Arow_ptr, Acol_ind, An, Brow_ptr, Bcol_ind, Bm, Csize)

Crow_ptr = int32(zeros(1,An+1));
Ccol_ind = int32(zeros(1,Csize));

coder.cinclude("BSpGEMM.h");
coder.cinclude("utils.h");

coder.ceval('SpGEMM_mat', coder.ref(Acol_ind), coder.ref(Arow_ptr), int32(An), ...
                          coder.ref(Bcol_ind), coder.ref(Brow_ptr), int32(Bm), ...
                          coder.ref(Ccol_ind), coder.ref(Crow_ptr));
           

end
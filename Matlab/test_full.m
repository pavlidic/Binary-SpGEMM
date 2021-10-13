function C = test_full(mat)
[A, ~, ~, ~] = mmread("./paralila 1h/"+mat+".mtx");
A = A';
tic;C=(A*A>0);toc
nnz(C)
end
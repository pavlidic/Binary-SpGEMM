function test_mtx(mat)
[A, ~, ~, ~] = mmread("./paralila 1h/"+mat+".mtx");
A=tril(A)';
tic;C=(A*A>0);toc
nnz(C)
end
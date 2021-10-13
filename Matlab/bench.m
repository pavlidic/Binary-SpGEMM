function bench()
n = 5e6;
d = 2; % approximate number of true elements per row

A = sprand( n, n, d/n ) > 0;
B = sprand( n, n, d/n ) > 0;

A=A*1;
B=B*1;

tic;
C = (A*B) > 0; 
toc
end
function write_spm(n,d,path)
%n = rows/collumns
%d = approximate number of true elements per row

A = sprand( n, n, d/n ) > 0;
%A = A*1;

mmwrite(path,A,'',"pattern");
end
function [Adash,err,CURerr]=CompressionTest(A,nsub,msub)
Adash=JacobiCompress(A,nsub,msub);
err=norm(A-Adash,'fro')/norm(A,'fro');
[C,U,R]=CUR(A,nsub,msub);
CURerr=norm(A-C*U*R,'fro')/norm(A,'fro');
disp(sprintf('JacobiError=%f CURerror=%f',err,CURerr));
end
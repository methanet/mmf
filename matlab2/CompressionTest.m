function [Adash,err,CURerr]=CompressionTest(A,nsub,msub,k)
if nargin<4 k=2; end
%Adash=JacobiCompress(A,nsub,msub);
Adash=MMFcompress(A,nsub,k);
err=norm(A-Adash,'fro')/norm(A,'fro');
[C,U,R]=CUR(A,nsub,msub);
CURerr=norm(A-C*U*R,'fro')/norm(A,'fro');
fprintf('JacobiError=%f CURerror=%f\n',err,CURerr);
end
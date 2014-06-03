function err=SymmCompressionTest(A,nrows,k)
if nargin<4 k=2; end
Adash=MMFcompress(A,nrows,k);
err=norm(A-Adash,'fro')/norm(A,'fro');
fprintf('MMFerror=%f ',err);
%[C,U,R]=CUR(A,nrows,nrows);
%CURerr=norm(A-C*U*R,'fro')/norm(A,'fro');
%fprintf('CURerror=%f ',CURerr);
[CN,W]=Nystrom(A,nrows);
NystrErr=norm(A-CN*W*CN','fro')/norm(A,'fro');
fprintf('NystromError=%f ',NystrErr);
fprintf('\n');
end
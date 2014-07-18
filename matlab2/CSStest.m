function [CSSerr, CSSrotErr] =CSStest(A,ncols,k)
if nargin<3, k=2; end
[C,X]=CSSrot(A,ncols,k);
CSSrotErr=norm(A-C*X,'fro')/norm(A,'fro');
fprintf('CSSrotError=%f ',CSSrotErr);
%[C,U,R]=CUR(A,nrows,nrows);
%CURerr=norm(A-C*U*R,'fro')/norm(A,'fro');
%fprintf('CURerror=%f ',CURerr);
[C2,X2]=CSS(A,ncols);
CSSerr=norm(A-C2*X2,'fro')/norm(A,'fro');
fprintf('CSSerror=%f ',CSSerr);
fprintf('\n');
end
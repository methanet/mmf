function [Adash,err]=CompletionTest(A,obs,Nrows,Ncols,k)
% [Adash,err]=CompletionTest(A,obs,nrows,ncols,k=2) 
% Recover a matrix from obs unifomrly sampled entries with MMF
% 
% A:           original matrix
% obs:         number of entries to sample
% nrows/ncols: number of rows/columns in MMF compressed form
% ncols:       order of rotations 
% Adash:       recovered matrix
% err:         relative Frobenius norm error
if nargin<5 k=2; end
[n,m]=size(A);
I=randi(n,obs,1);
J=randi(m,obs,1);
v=zeros(obs,1);
for o=1:obs
    v(o)=A(I(o),J(o));
end
[U,S,V]=MMFcomplete(I,J,v,n,m,Nrows,Ncols,k);
Adash=U*S*V';
%A-Adash
err=norm(A-Adash,'fro')/norm(A,'fro');
fprintf('JacobiError=%f\n',err);
end
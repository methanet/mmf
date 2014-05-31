function [Adash,err, errOmega]=SymMMFCompletionTest(A,obs,Nrows,Ncols,k)
% COMPLETIONTEST attempt to recover a matrix from sampled entries with MMF
% 
% Usage: [Adash,err]=CompletionTest(A,obs,nrows,ncols,k=2) 
%
% A:           orginial matrix
% obs:         number of entries to sample
% nrows/ncols: number of rows/columns in MMF compressed form
% ncols:       order of rotations 
% Adash:       recovered matrix
% err:         relative Frobenius norm error
%
if nargin<5 k=2; end
[n,m]=size(A);
%I= sparse(randi(n,obs,1));
%J= sparse(randi(m,obs,1));

%added
[nnzrows, nnzcols] = find(A); 
[nnzsize,~] = size(nnzrows);
omega = randsample(nnzsize, obs); 
I= nnzrows(omega);
J = nnzcols(omega);


v=zeros(obs,1);
for o=1:obs
    v(o)=A(I(o),J(o));
end

[U,S,V]=SymMMFcomplete(I,J,v,n,Nrows,k);
Adash=U*S*V';

err=norm(A-Adash,'fro')/norm(A,'fro');
errOmega = norm(A (sub2ind(size(A), I, J))-Adash(sub2ind(size(A), I, J)), 'fro')/norm(A(sub2ind(size(A), I, J)) ,'fro'); 

fprintf('SymMMFJacobiError=%f\n',err);
end

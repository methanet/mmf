function [U,S,V]=MMFcomplete(I,J,v,n,m,nrows,ncols,k)
% [U,S,V]=MMFcomplete(I,J,v,n,m,Nrows,Ncols,k) 
% Multiresolution matrix completion
%
% I,J:   vector of row resp. column indices of observed entries
% v:     vector of observed entries
% n,m:   dimesions of full matrix
% nrows: number of rows of compressed matrix
% ncols: number of columns of compressed matrix
% k=2:   order of rotations
% U,S,V: estimated matrix is returned in the form U*S*V'
if nargin<8 k=2; end
preA=zeros(n,m);
[obs,~]=size(I);
for o=1:obs 
    preA(I(o),J(o))=v(o);
end;
if k==2 
    [U,A,rows]=RowJacobi(preA,nrows);
    [V,A,cols]=RowJacobi(preA',ncols);
else
    [U,A,rows]=RowMMF(preA,nrows,k);
    [V,A,cols]=RowMMF(preA',ncols,k);
end
U=U(:,rows);
V=V(:,cols);
M=zeros(obs,nrows*ncols);
for o=1:obs
    for i=1:nrows
        for j=1:ncols
            p=(i-1)*ncols+j;
            M(o,p)=U(I(o),i)*V(J(o),j);
        end
    end
end
x=pinv(M)*v;
%F=M*x;
%F(1:10)
%v(1:10)
%[mu,ms,mv]=svd(M);
%diag(ms)
S=reshape(x,nrows,ncols)';
end
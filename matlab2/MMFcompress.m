function Adash=MMFcompress(A,nrows,k)
% Adash=MMFcompress[A,nrows,k]
% Compress a symmetric matrix to an nrows*nrows core plus diagonals with MMF
% 
% A:     input matrix
% nrows: size of core
% k:     order of rotations
[n,~]=size(A);
[U,A,active]=MMF(A,nrows,k);
Adash=zeros(n,n);
for i=1:n Adash(i,i)=A(i,i); end
S=A(active,active);
Adash(active,active)=S;
Adash=U*Adash*U';

function [U,A,active]=RowMMF(A,nrows,k)
% [U,Adash,active]=RowMMF(A,nrows,k):   Compute the one-sided (row-wise) 
% Jacobi Multiresolution Matrix Factorization of the symmetric matrix A
%
% A:      the input matrix
% nrows:  dimensionality of final core matrix
% k:      order of final core matrix
% U:      orthogonal matrix, each column of which is a wavelet
% Adash:  transformed version of A
% active: list of rows of A still active at the end of the process
%
[n,~]=size(A);
active=1:n;
U=eye(n);

for r=n:-1:nrows+1
    ri=randi(r);
    i=active(ri);
    inp=A(active,:)*A(i,:)';
    inp(ri)=-Inf;
    [~,Ir]=kmax(inp,k);
    I=active(Ir);
    S=A(I,:)*A(I,:)';
    [V,D]=eig(S);
    U(:,I)=U(:,I)*V;
    A(I,:)=V'*A(I,:);
    [~,i]=min(diag(D));
    active=[active(1:Ir(i)-1);active(Ir(i)+1:r)];
    fprintf('%d %f\n',I(i),D(i,i));
end
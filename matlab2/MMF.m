function [U,A,active]=MMF(A,nrows,k)
% [U,Adash,active]=MMF(A,nrows,k):   Compute the Jacobi Multiresolution 
% Matrix Factorization of the symmetric matrix A
%
% A:      the input matrix
% nrows:  dimensionality of final core matrix
% k:      order of final core matrix
% U:      orthogonal matrix, each column of which is a wavelet
% Adash:  transoformed version of A
% active: list of rows/columns of A that form the final core

[n,~]=size(A);
active=1:n;
U=eye(n);

for r=n:-1:nrows+1
    ri=randi(r);
    i=active(ri);
    inp=A(active,active)*A(i,active)';
    inp(ri)=-Inf;
    [~,Ir]=kmax(inp,k);
    I=active(Ir);
    S=A(I,active)*A(I,active)';
    [V,D]=eig(S);
    U(:,I)=U(:,I)*V;
    A(I,:)=V'*A(I,:); 
    A(:,I)=A(:,I)*V;
    [~,i]=min(diag(D));
    active=[active(1:Ir(i)-1),active(Ir(i)+1:r)];
    for j=1:k fprintf('%d ',I(j)); end 
    fprintf('%f\n',2*D(i,i));
    %fprintf('%d %f\n',I(i),2*D(i,i));
end

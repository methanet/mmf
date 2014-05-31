function [U,S,V]=MMFcomplete(I,J,v,n,m,Nrows,Ncols,k)
if nargin<8 k=2; end
preA=zeros(n,m);
[obs,~]=size(I);
for o=1:obs 
    preA(I(o),J(o))=v(o);
end;
if k==2 
    [U,A,rows]=RowJacobi(preA,Nrows);
    [V,A,cols]=RowJacobi(preA',Ncols);
else
    [U,A,rows]=RowMMF(preA,Nrows,k);
    [V,A,cols]=RowMMF(preA',Ncols,k);
end
U=U(:,rows);
V=V(:,cols);
M= zeros(obs,Nrows*Ncols);

for o=1:obs
    for i=1:Nrows
        for j=1:Ncols
            p=(i-1)*Ncols+j;
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

S=sparse(reshape(x,Nrows,Ncols)');

end
function [U,S,V]=SymMMFcomplete(I,J,v,n,Nsubs, k)
if nargin<6 k=2; end
preA=zeros(n,n);
[obs,~]=size(I);
for o=1:obs 
    preA(I(o),J(o))=v(o);
end;
preA*preA'
preA'*preA
if k==2 
    [U,~,rows]=JacobiMMF(preA*preA',Nsubs);
    [~,V,cols]=JacobiMMF(preA'*preA,Nsubs);
else
    [U,~,rows]=MMF(preA*preA',Nsubs,k);
    [~,V,cols]=MMF(preA'*preA,Nsubs,k);
end
rows
cols
U
V
U=U(:,rows);
V=V(:,cols);

M= zeros(obs,Nsubs*Nsubs);

for o=1:obs
    for i=1:Nsubs
        for j=1:Nsubs
            p=(i-1)*Nsubs+j;
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

S=sparse(reshape(x,Nsubs,Nsubs)');

end

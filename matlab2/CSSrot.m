function [C,X]=CSSrot(A,ncols,k)
[~,m]=size(A);
S=A'*A;
S=abs(S-diag(diag(S)));
knn=zeros(m,k);
p=zeros(m,1);
for i=1:m
    [~,I]=kmax(S(i,:),k); 
    knn(i,:)=I;
    p(i)=min(eig(S(I,I)));
end
counts=mnrnd(ncols,p/sum(p));
selector=[];
for i=1:m
    selector=[selector,i*ones(1,counts(i))];
end
C=A(:,selector);
X=pinv(C'*C)*C'*A;


    

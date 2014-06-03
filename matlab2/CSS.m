function [C,X]=CSS(A,ncols)
[~,m]=size(A);
p=sum(A.^2,1);
selector=[];
counts=mnrnd(ncols,p/sum(p));
selector=[];
for i=1:m
    selector=[selector,i*ones(1,counts(i))];
end
C=A(:,selector);
X=pinv(C'*C)*C'*A;
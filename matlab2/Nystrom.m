function [C,W]=Nystrom(A,nrows)
[n,~]=size(A);
selector=randi(n,nrows,1);
C=A(:,selector);
W=pinv(A(selector,selector));
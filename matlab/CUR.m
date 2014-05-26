function [C,U,R]=CUR(A,nsub,msub)
[n,m]=size(A);
Cselector=randi(m,msub,1);
C=A(:,Cselector);
Rselector=randi(n,nsub,1);
R=A(Rselector,:);
U=(n/nsub)*pinv(C'*C)*C(Rselector,:)';
end
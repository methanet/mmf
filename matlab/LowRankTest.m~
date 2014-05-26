function [A,Adash,err,CURerr]=LowRankTest(n,m,k,nsub,msub)
U=2*rand(n,k)-1.0;
N=sqrt(sum(U.^2)).^(-1);
U=U.*(ones(n,1)*N);
V=2*rand(m,k)-1.0;
N=sqrt(sum(V.^2)).^(-1);
V=V.*(ones(m,1)*N);
A=U*V';
Adash=JacobiCompress(A,nsub,msub);
err=norm(A-Adash,'fro')/norm(A,'fro');
[C,U,R]=CUR(A,nsub,msub);
CURerr=norm(A-C*U*R,'fro')/norm(A,'fro');
disp(sprintf('JacobiError=%f CURerror=%f',err,CURerr));
end
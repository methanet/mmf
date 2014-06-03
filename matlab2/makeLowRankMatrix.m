function A=makeLowRankMatrix(n,m,k)
U=2*rand(n,k)-1.0;
N=sqrt(sum(U.^2)).^(-1);
U=U.*(ones(n,1)*N);
V=2*rand(m,k)-1.0;
N=sqrt(sum(V.^2)).^(-1);
V=V.*(ones(m,1)*N);
A=U*V';
end
function [U,A,activelist]=RowMMF(A,Nrows,k)
[n,~]=size(A);
activelist=zeros(n,1); for i=1:n activelist(i)=i; end
U=eye(n);

for r=n:-1:Nrows+1
    ri=randi(r);
    i=activelist(ri);
    inp=A(activelist,:)*A(i,:)';
    inp(ri)=-Inf;
    [~,Ir]=kmax(inp,k);
    I=activelist(Ir);
    S=A(I,:)*A(I,:)';
    [V,D]=eig(S);
    U(:,I)=U(:,I)*V;
    A(I,:)=V'*A(I,:);
    [~,i]=min(diag(D));
    activelist=[activelist(1:Ir(i)-1);activelist(Ir(i)+1:r)];
    fprintf('%d %f\n',I(i),D(i,i));
end
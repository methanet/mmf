function MRF(A,nsub)
[n,m]=size(A);

T=zeros(n);
for i=1:n T(i)=norm(A(i,:))^2; end

S=zeros(n,n);
for i=1:n
    for j=1:i-1
        Tij=A(i,:)*A(j,:)';
        S(i,j)=T(i)+T(j)-sqrt((T(i)+T(j))^2-4*T(i)*T(j)+4*Tij^2);
        S(j,i)=S(i,j);
    end
end

activelist=zeros(n);
for i=1:n activelist(i)=i; end

U=eye(n);
for r=n:-1:sub
    ri=rand(r);
    i=activelist(ri);
    [dummy,j]=min(S(i,:));
    tau=(T(j)-T(i))/(A(i,:)*A(j,:)')/2;
    t=-tau+sqrt(tau^2+1);
    c=1/sqrt(1+t^2);
    Asub=A([i,j],:);
    R=[[c,c*t];[-c*t,c]];
    Ar=R*Asub;
    A(i,:)=Ar(1,:);
    A(j,:)=Ar(2,:);
    U()
    activelist=[activelist(1:ri-1),activelist(ri+1,r)];
    
end


end
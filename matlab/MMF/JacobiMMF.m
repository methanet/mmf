function [U,A,activelist]= JacobiMMF(A,nsub)
[n,m]=size(A);

T=zeros(n);
for i=1:n T(i)=norm(A(i,:))^2; end

S=zeros(n,n);
for i=1:n
    for j=1:i-1
        Tij=A(i,:)*A(j,:)';
        S(i,j)=(T(i)+T(j)-sqrt((T(i)+T(j))^2-4*T(i)*T(j)+4*Tij^2))/2;
        S(j,i)=S(i,j);
    end
end

activelist=zeros(n,1);
for i=1:n activelist(i)=i; end

U=eye(n);
scores=zeros(n,1);
for r=n:-1:nsub+1 
    ri=randi(r);
    i=activelist(ri);
    r
    for rj=1:r
        j=activelist(rj);
	score(rj)=(T(i)+T(j)-sqrt((T(i)+T(j))^2-4*T(i)*T(j)+4*(A(i,:)*A(j,:)')^2));
	(T(i)+T(j)-sqrt((T(i)+T(j))^2-4*T(i)*T(j)+4*(A(i,:)*A(j,:)')^2))/2
    end
    ri
    score
    score(ri)=NaN;
    score
    [dummy,rj]=min(score(1:r));
    
    %[dummy,rj]=min(S(i,:));

    j=activelist(rj);
    tau=(T(j)-T(i))/(A(i,:)*A(j,:)')/2;
    t=1/(abs(tau)+sqrt(tau^2+1));
    c=1/sqrt(1+t^2);
    Asub=A([i,j],[i,j]);
    R=[[c,-c*t];[c*t,c]];

    A([i,j],:)=R*A([i,j],:);
    A(:,[i,j])=A(:,[i,j])*R';
    U(:,[i,j])=U(:,[i,j])*R';

   
    if (norm(A(j,:))>norm(A(i,:)))
	T(j)=T(i)+T(j)-score(rj);
        activelist=[activelist(1:ri-1);activelist(ri+1:r)];
    else 
	T(i)=T(i)+T(j)-score(rj);
        activelist=[activelist(1:rj-1);activelist(rj+1:r)];
    end
    
end
U
A
end

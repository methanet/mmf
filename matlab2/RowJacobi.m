function [U,A,activelist]=RowJacobi(A,nsub)
[n,m]=size(A);
T=zeros(n); for i=1:n T(i)=norm(A(i,:))^2; end
activelist=zeros(n,1); for i=1:n activelist(i)=i; end

U=eye(n);
scores=zeros(n,1);
for r=n:-1:nsub+1
    ri=randi(r);
    i=activelist(ri);
    for rj=1:r
        j=activelist(rj);
        score(rj)=(T(i)+T(j)-sqrt((T(i)+T(j))^2-4*T(i)*T(j)+4*(A(i,:)*A(j,:)')^2))/2;
    end
    score(ri)=NaN;
    [dummy,rj]=min(score(1:r));
    j=activelist(rj);
    %score(rj)
    %[T(i),A(i,:)*A(i,:)'] % OK
    %eig(A([i,j],:)*A([i,j],:)')
    prenm=[norm(A(i,:)),norm(A(j,:))];
    tau=(T(j)-T(i))/(A(i,:)*A(j,:)')/2;
    %tau=2*A(i,:)*A(j,:)'/(T(j)-T(i));
    t=1/(abs(tau)+sqrt(tau^2+1));
    if tau<0 t=-t; end 
    c=1/sqrt(1+t^2);
    R=[[c,-c*t];[c*t,c]];
    A([i,j],:)=R*A([i,j],:); 
    U(:,[i,j])=U(:,[i,j])*R';
    if(norm(A(j,:))>norm(A(i,:)))
        T(j)=T(i)+T(j)-score(rj);
        activelist=[activelist(1:ri-1);activelist(ri+1:r)];
        postnm=[norm(A(i,:)),norm(A(j,:))];
        Tj=T(j);
        display(sprintf('%d %d %f',j,i,score(rj)));
    else
        T(i)=T(i)+T(j)-score(rj);
        activelist=[activelist(1:rj-1);activelist(rj+1:r)];
        postnm=[norm(A(j,:)),norm(A(i,:))];
        Tj=T(i);
        display(sprintf('%d %d %f',i,j,score(rj)));
    end
    %display(sprintf('%f %f %f %f %f',prenm(1)^2,prenm(2)^2,postnm(1)^2,postnm(2)^2,Tj));
    % display(sprintf('%f %f',norm(prenm),norm(postnm))); % OK
end


end
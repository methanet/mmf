function A=makeNetworkMatrix(E)
[e,~]=size(E);
maxlabel=max(max(E(:,1)),max(E(:,2)));
flag=zeros(maxlabel,1);
flag(E(:,1))=1;
flag(E(:,2))=1;
t=1;
dict=zeros(maxlabel,1);
for i=1:maxlabel
    if flag(i)==1 
        dict(i)=t;
        t=t+1;
    end
end
n=t-1;
A=zeros(n,n);
for i=1:e
    A(dict(E(i,1)),dict(E(i,2)))=1;
    A(dict(E(i,2)),dict(E(i,1)))=1;
end
    
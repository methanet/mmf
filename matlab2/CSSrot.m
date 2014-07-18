function [C,X]=CSSrot(A,ncols,k)
[~,m]=size(A);
S=A'*A;
S=abs(S-diag(diag(S)));
knn=zeros(m,k);
p=zeros(m,1);
%selector=[];
globalselector = []; 

iterCount = 0;
while numel(unique(globalselector)) < m-ncols
    iterCount =0;
    for i=1:m
        [~,I]=kmax(S(i,:),k);
        knn(i,:)=I;
        p(i)= 1/min(eig(S(I,I)));
    end
    counts=mnrnd(m, p/sum(p));
    indexselector=[];
    for j=1:m
        %selector=[selector,i*ones(1,counts(j))];
        [~,Inds] = sort(counts, 'descend');
        indexselector = Inds(1:200);
    end
    
    for ind = 1:numel(indexselector)
        cols = knn (indexselector(ind),:);
        SS = A(:,cols)'*A(:,cols); 
        [V,~]=eig(SS);
        A(:, cols)= A(:, cols)*V;
    end 
    globalselector = [globalselector, indexselector];
end
size(globalselector)
active = setdiff(1:m, globalselector);
size(active)
C=A(:, active);
X=pinv(C'*C)*C'*A;
CSSrotErr=norm(A-C*X,'fro')/norm(A,'fro');
CSSrotErr
fprintf ('Iterations: %d', iterCount);


    

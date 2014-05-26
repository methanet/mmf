function Adash=JacobiCompress(A,nsub,msub)
[U,A1,alist1]=RowJacobi(A,nsub);
[V,A2,alist2]=RowJacobi(A1(alist1,:)',msub);
Adash=U(:,alist1)*(V(:,alist2)*A2(alist2,:))';
end
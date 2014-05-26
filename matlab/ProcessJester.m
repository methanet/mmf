function A=ProcessJester(A)
[I,J]=find(A==99);
[n,dummy]=size(I);
for i=1:n
    A(I(i),J(i))=0;
end
end
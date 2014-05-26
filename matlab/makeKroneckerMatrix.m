function A=makeKroneckerMatrix(G,k)
A=G;
for i=1:k-1
    A=kron(A,G);
end
end
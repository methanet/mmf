%obs is fraction e.g. 0.7
function phasetransitionKronecker(obs,iters, step, kronsize)%, Nrows, Ncols)
generator = rand(2,2); 
generator = triu(generator)' + triu(generator,1); 
M = makeKroneckerMatrix (generator, kronsize);
imagesc(M);
[n1,n2] = size(M);
size(M)
%oversampling = 5; 
m = floor(n1*n2*obs);%min(5*df,round(.99*n1*n2) ); i

MRFerrorOmega = zeros(ceil(n1/step),1);
MRFerror = zeros(ceil(n1/step),1); 

for i=1: iters
    count =1;
    for j= 1:step:n1
        A = M; 
        [~,err, errOmega]=CompletionTest(A,floor(m),j,j); 
        MRFerrorOmega (count) = MRFerrorOmega (count)+ errOmega;
        MRFerror (count) =  MRFerror (count)+ err; 
        count= count+1;
    end
end

plot(1:step:n1, MRFerror/iters, 'b', 1:step:n1, MRFerrorOmega/iters, 'r')
legend('MRFerror', 'MRFerrorOmega') ; 
end
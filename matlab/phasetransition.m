function phasetransition(obs,iters, kronsize)%, Nrows, Ncols)
probs = 0.1:0.05:0.9; 
[~,n] = size(probs);
MRFerrorOmega = zeros(n,1);
MRFerror = zeros(n,1);

for i=1: iters
    parfor  p = 1:n
        generator = [[1,probs(p)]; [probs(p), 1]];
        generator = triu(generator)' + triu(generator,1); 
        M = makeKroneckerMatrix (generator, kronsize);
        [n1,n2] = size(M);
        m = floor(n1*n2*obs);%min(5*df,round(.99*n1*n2) ); i
        [~,err, errOmega]=CompletionTest(M,m,n1/2,n1/2); 
        MRFerrorOmega (p) = MRFerrorOmega (p)+ errOmega;
        MRFerror (p) =  MRFerror (p)+ err;
    end
end
plot(probs, MRFerror/iters, 'b', probs, MRFerrorOmega/iters, 'r')
title(sprintf('Kronecker 2x2 probabilities vs error, size: %d',2^kronsize)); 
legend('MRFerror', 'MRFerrorOmega') ; 
xlabel('probability[[1,p];[p,1]]')
end
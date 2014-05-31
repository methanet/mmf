%obs is the fraction  (i.e 0.5, 0.7 etc) of entries to sample
% kronsize determines the size of the final Kronecker graph 
%SVT params: rank r, 
%MMF params: nrows, ncols, 
function compareKronecker(obs, iters, kronsize, r, Nrows, Ncols)
generator = rand(2,2); 
generator = triu(generator)' + triu(generator,1); 
M = makeKroneckerMatrix (generator, kronsize);
imagesc(M);
[n1,n2] = size(M);
size(M)
df = r*(n1+n2-r);
%oversampling = 5; 
m = floor(n1*n2*obs);%min(5*df,round(.99*n1*n2) ); i

SVTerrorOmega = zeros(iters,1);
SVTerror = zeros(iters,1); 

MRFerrorOmega = zeros(iters,1);
MRFerror = zeros(iters,1); 

for i=1: iters
    Omega = randsample(n1*n2, m); %observed entries
    % nnzero = find(M==1); 
    % obs = 0.7*numel(nnzero); 
    % Omega = randperm(numel(M)); Omega = Omega(1:obs);

    data = M(Omega);
    %sigma = 0.0 ; 
    %sigma = .05*std(data);
    %data = data + sigma*randn(size(data));
    tau = 5*sqrt(n1*n2); 
    delta = 1.5; 
    fprintf('Percent sampled entries: %d\n', floor(m*100/(n1*n2)));

    maxiter = 500; 
    tol = 1e-3;
    [U,S,V,numiter] = SVT([n1 n2],Omega,data,tau,delta,maxiter,tol);    
    X = U*S*V';
    fprintf('The recovered rank is %d\n',length(diag(S)) );
    fprintf('The relative error on Omega is: %d\n', norm(data-X(Omega))/norm(data))
    fprintf('The relative recovery error is: %d\n', norm(M-X,'fro')/norm(M,'fro'))
    fprintf('The relative recovery in the spectral norm is: %d\n', norm(M-X)/norm(M))
    fprintf('Actual  number of iterations %d\n',numiter);
    SVTerrorOmega (i) =  norm(data-X(Omega))/norm(data); 
    SVTerror (i) = norm(M-X,'fro')/norm(M,'fro');  
end

for i=1: iters
    A = M; 
    imagesc(A)
    [~,err, errOmega]=CompletionTest(A,floor(m),Nrows,Ncols); 
    MRFerrorOmega (i) =  errOmega;
    MRFerror (i) = err;
  
end

fprintf('SVT: error on Omega %d: error total: %d\n',mean (SVTerrorOmega), mean (SVTerror))
fprintf('MMF: error on Omega %d: error total: %d\n',mean (MRFerrorOmega), mean (MRFerror))


function CSSgraph(A,iters, filename, k)
if nargin<4 k=2; end
cols = 3:2:25;
CSSrotError = zeros (iters,numel(cols));
CSSerror = zeros (iters, numel(cols));

for i=1:iters
	for j=1:numel(cols)
		[CSSerr, CSSrotErr]=CSStest(A, cols(j), k);
		CSSerror(i,j) = CSSerr;
		CSSrotError(i,j) = CSSrotErr; 
	end
end
meanCSSerror = mean(CSSerror,1); 
meanCSSrotError = mean(CSSrotError,1);
stdCSSerror = std(CSSerror,0,1); 
stdCSSrotError = std(CSSrotError,0,1);

f=figure;
%errorbar( cols, meanCSSerror, stdCSSerror, 'b');
%hold on;
%errorbar(cols, meanCSSrotError,stdCSSrotError, 'r');
plot (cols, meanCSSerror,'b', cols, meanCSSrotError, 'r'); 
title(sprintf('iters=%d, k=%d',iters, k));
legend('CSSerror', 'CSSrotError') ;
xlabel('columns'); 
print(f,'-djpeg', filename);
end	

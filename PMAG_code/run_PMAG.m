function [results,alpha,runtime,obj] = run_PMAG(X,Y,projDim,r)
[Y, idx] = sort(Y);
X = X(idx, :);
X = full(double(X));  % n*d
if min(Y)==0
    Y = Y + 1;
end
c = length(unique(Y));
countC = histcounts(Y);
fprintf('The clusters number are: %d \n',countC);

%%
options.PCARatio = 0.9;
[eigvector, ~] = PCA(X,options); %   PCA  Ratio
X = X*eigvector;
[W,alpha,obj,runtime] = PMAG_Eproj(X,projDim,c,r);
F = kmeans(real(X*W),c);
results = Clustering8Measure(Y, F);
end


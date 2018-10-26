

function [Xnew, ind] = cleanDesignMatrix(X);

% clean up the design matrix by removing constant columns.
v  = var(X, [], 1);
ind = find(v==0);
Xnew = X( :, setdiff(1:size(X,2), ind));


function [maxCC, maxLag, Consistent] = computeSimilarityPeaks( X, Y )
[cc,lag] = xcorr( X, Y, 'coeff');
[maxCC, ind] = max(cc);
 maxLag = lag(ind);

if maxLag < 5 && maxLag > -5
    Consistent = 1;
else
    Consistent = 0;
end;
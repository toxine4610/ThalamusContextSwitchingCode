
function [ExpH, perf1, Toss] = getBehavior(Z_context1)


Toss = Z_context1(:,1);

w = gausswin(5);
w = w/sum(w);
perf1 =conv(Z_context1(:,1),w,'same');
[ExpH, rup, rlo ] = computeExpectation(Z_context1);
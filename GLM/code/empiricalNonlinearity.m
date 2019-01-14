function [sprate, binCenters, praw]=empiricalNonlinearity(st,xproj,nBins) %% check nonlinearity
% [pspk, binCenters]=empiricalNonlinearity(st,xproj,nBins) %% check nonlinearity

if ~exist('nBins', 'var')
    nBins    = 50;
end

binEdges   = (quantile(xproj, linspace(0,1,nBins+1)));
% binCenters = (quantile(xproj, linspace(1/(2*nBins),1-(1/(2*nBins)), nBins)));

binEdges(diff(binEdges)<=0)=[];
binCenters=binEdges(1:end-1)+diff(binEdges)/2;
% binCenters(diff(binCenters)<=0)=[];
% histogram xproj and kTxspk
[praw, id] = histc(xproj, binEdges); 
praw = praw(1:numel(binCenters)); praw(end) = [];
binCenters(praw==0) = [];

ids = unique(id);
sprate = zeros(numel(ids)-1,1);

for k = 1:numel(ids)-1
    sprate(k) = mean(st(id==ids(k)));
end

nb=numel(binCenters);
ns=numel(sprate);
np=numel(praw);
if (nb ~= ns) || (nb ~= np)
    mn=min([nb ns np]);
    binCenters=binCenters(1:mn);
    sprate=sprate(1:mn);
    praw=praw(1:mn);
end
    
% pspk = sprate./praw;
function [bc,mm,sd,nn,lm]=errorbarBin(xx,yy,bins, varargin)
% bin x errorbar y
if ~exist('bins', 'var') || isempty(bins)
    bins=linspace(min(xx), max(xx), 10);
end

mm=arrayfun(@(x,y) nanmean(yy(xx>x & xx<y)), bins(1:end-1), bins(2:end));
sd=arrayfun(@(x,y) nanstd(yy(xx>x & xx<y)), bins(1:end-1), bins(2:end));
nn=arrayfun(@(x,y) nansum(xx>x & xx<y), bins(1:end-1), bins(2:end));

lm=fitlm(xx,yy);
sd=sd./sqrt(nn);
bc=bins(1:end-1)+diff(bins)/2;
errorbar(bc, mm, sd, varargin{:})
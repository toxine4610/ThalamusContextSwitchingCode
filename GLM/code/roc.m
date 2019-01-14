function [AUC] = roc(Xhat, Xtrue, alpha)
% [AUC] = roc(Xhat, Xtrue)
% calculate the area under the ROC curve for Xhat with labels Xtrue
if nargin < 3
    alpha = .05; 
    if nargin < 2
        help roc
    end
end

Xhat = Xhat(:);
Xtrue = Xtrue(:);

% xlevels = linspace(min(Xhat), max(Xhat), 20);
xlevels = [-inf; unique(Xhat)-mean(diff(Xhat)); inf];
nlevels = numel(xlevels);
yroc = zeros(nlevels,1);
xroc = zeros(nlevels,1);
nT = sum(Xtrue);
nF = sum(~Xtrue);

for ii = 1:nlevels   
   tot = Xhat>xlevels(ii);
   yroc(ii) = sum(Xtrue & tot)/nT;
   xroc(ii) = sum(~Xtrue & tot)/nF;
end

% remove excess zeros -- these can mess up the trapz method of finding area
% under a curve
xz = find(xroc==0);
[my, yi] = max(yroc(xz));

support = [find(xroc ~=0); xz(yi)];
% plot(xroc(support), yroc(support), 'o-r', [0 1], [0 1], 'k')

Area=abs(trapz(xroc(support),yroc(support))); %estimate the area under the curve
%standard error of area
Area2=Area^2; Q1=Area/(2-Area); Q2=2*Area2/(1+Area);
V=(Area*(1-Area)+(nF-1)*(Q1-Area2)+(nT-1)*(Q2-Area2))/(nF*nT);
Serror=realsqrt(V);
%confidence interval
ci=Area+[-1 1].*(realsqrt(2)*erfcinv(alpha)*Serror);
if ci(1)<0; ci(1)=0; end
if ci(2)>1; ci(2)=1; end
%z-test
SAUC=(Area-0.5)/Serror; %standardized area
p=1-0.5*erfc(-SAUC/realsqrt(2)); %p-value

AUC.AUC = Area;
AUC.serror = Serror;
AUC.ci = ci;
AUC.p = p;
AUC.xr = xroc; 
AUC.yr = yroc;
function [T1,T2, expl1, expl2] = extractFilterPC( CouplingFilter, Lc )


%%
warning off;

clear x
LL = Lc./1e4;
ind = find(LL>0.3);
ct  = 0;

for i = 1:length(ind)
    ct=ct+1;
    x(ct,:,:) = cell2mat( CouplingFilter{ind(i)} )';
end;

x = permute( x, [2,1,3] ); % numMD x numPFC x time


for i = 1:size(x,1)
    foo = squeeze( x(i,:,:) );

    [coeff,score,latent,tsquared,explained,mu] = pca(foo);
    T1(i,:) = exp(coeff(:,1))';
    T2(i,:) = exp(coeff(:,2))';
    
    expl1(i) = explained(1);
    expl2(i) = explained(2);
end;



%%% Computes the principal components of MD and PFC inputs into PFC cells.

function [V,MDInput1, MDInput2, PFCInput1, PFCInput2, SNRMD, SNRPFC, MD_expl1, MD_expl2, PFC_expl1, PFC_expl2] = classifyInputs(Lc, CouplingFilter, CouplingFilterPFC)

%%
warning off;

clear x
% LL = Lc./1e4;
% ind = find(LL>0.3);

ind = 1:size(CouplingFilterPFC,2);

%%

ct  = 0;

for i = 1:length(ind)
    ct=ct+1;
    clear x
    %==== PCA of MD inputs to this cell.
    x  = cell2mat( CouplingFilterPFC{ind(i)} )';
    V(i,i) = NaN;
    indexes = setdiff( 1:(size(x,1)+1), i);
    for j = 1:length(indexes)
        V(i,indexes(j)) =  trapz( x(j,1:20) );
    end;
end;
    

%%

ct  = 0;
for i = 1:length(ind)
    ct=ct+1;
    clear x
    %==== PCA of MD inputs to this cell.
    x  = cell2mat( CouplingFilter{ind(i)} )';
    clear coeff score explained
    [coeff,score,latent,tsquared,explained,mu] = pca(x);
    MDInput1(ct,:) = exp(coeff(:,1))';
    MDInput2(ct,:) = exp(coeff(:,2))';
    
    MD_expl1(ct) = explained(1);
    MD_expl2(ct) = explained(2);
    
    SNRMD(ct) = mean( MDInput1(ct,1:10) )./mean( MDInput1(ct,50:70) );

    
    clear x
    %==== PCA of PFC inputs to this cell.
    x  = cell2mat( CouplingFilterPFC{ind(i)} )';
    clear coeff score explained
    [coeff,score,latent,tsquared,explained,mu] = pca(x);
    PFCInput1(ct,:) = exp(coeff(:,1))';
    PFCInput2(ct,:) = exp(coeff(:,2))';
    
    PFC_expl1(ct) = explained(1);
    PFC_expl2(ct) = explained(2);
    
    SNRPFC(ct) = mean( PFCInput1(ct,1:10) )./mean( PFCInput1(ct,50:70) );
    
end;

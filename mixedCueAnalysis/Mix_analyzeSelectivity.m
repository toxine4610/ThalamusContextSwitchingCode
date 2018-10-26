


foo = cell2mat( PFC_dp_C1' );
foo = foo(:, indDelay );
muDP1 = nanmean(foo,2);


foo2 = cell2mat( PFC_dp_C2' );
foo2 = foo2(:, indDelay );
muDP2 = nanmean(foo2,2);

ind1 = find( isnan(muDP1)==0);
ind2 = find( isnan(muDP2)==0);
indkeep = intersect(ind1, ind2);


p = polyfit( muDP1(indkeep), muDP2(indkeep), 1);
mdl = polyval(p, muDP1);

%%

TS  = foo( indkeep,:);
TS2 = foo2(indkeep,:);


[coeff,score,latent,tsquared,explained,mu] = pca(TS);
[coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(TS2);

for i = 1:12
    this = coeff(:, i);
   [v, ind] = max( this );
   t2p(i)   = ind;
end

[~, sorting_c1] =  sort( t2p );

clear t2p 
for i = 1:12
    this = coeff2(:, i);
   [v, ind] = max( this );
   t2p(i)   = ind;
end

[~, sorting_c2] =  sort( t2p );



figure(1);
cmap = hot(12);
for i = 1:12
   plot( time(indDelay), coeff(:,sorting_c1(i))','color',cmap(i,:),'linewidth',3); hold on
end

figure(2);
% cmap = hot(12);
for i = 1:12
   plot( time(indDelay), coeff2(:,sorting_c2(i))','color',cmap(i,:),'linewidth',3); hold on
end


%%

for i = 1:12
    sc1 = score(:, sorting_c1(i) );
    sc2 = score(:, sorting_c2(i) );
    indpos1 = intersect( find(sc1>0.5), find(sc2>0.5) );
    figure(i);
    scatter( sc1(indpos1), sc2(indpos1));
    IND{i} = indpos1;
end;


function compBest = svdDecompWaveforms(waveforms_this_unit)

rng(1);

%%
for i = 1:size(waveforms_this_unit,2)
    foo = squeeze( waveforms_this_unit(:, i, :));
    [u1,s1,v1] = svd(double(foo),'econ');
    firstcomp(i,:) = u1(:,1)';
    sv(i)          = s1(1);
    
    r(i)  = range(firstcomp(i,:));
end;
[svsort, ind] = sort(sv);
[~, bestind]  = max(sv);

ind = find( r ~= 1);
%%

clear CC

P = nchoosek( ind ,2);
for i=1:size(P,1)
  CC(i) = corr( firstcomp( P(i,1), :)', firstcomp( P(i,2), :)' );
end;

[vvv, iii] = sort(CC);
P1 = P(iii(end), 1); P2 = P(iii(end), 2);

compBest      = [firstcomp(P1, :); firstcomp(P2,:)];
compBest      = mean(compBest,1) ;
compBest       = sgolayfilt(compBest, 3, 7);
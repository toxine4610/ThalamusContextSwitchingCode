
ct = 0;
for i = 2:length(ISI)-1
    
    ct = ct+1;
    currSpike = spiketimes_thiscell(i);
    beforeSpike = spiketimes_thiscell(i-1);
    afterSpike = spiketimes_thiscell(i+1);
    
    ISIbefore(ct) = currSpike-beforeSpike;
    ISIafter(ct)  = afterSpike-currSpike;
    
end;

%%

X = [(ISIafter);(ISIbefore)];
% indexkeep = intersect( find(X(1,:)<-2.5), find( X(2,:)< -2.5) );
% X = X(:, indexkeep);

rng(1);
opts = statset('Display','final');
[idx,C] = kmeans(X',2,'Distance','cityblock','Replicates',25,'Options',opts);

function [spikeRate, spikeMat, timeOut, indBase] = makeSpikeRates(foo, range , bin, filtWidth)

clCell = foo;

time = range(1):bin:range(2);
filtSize = round(filtWidth./bin);
filter = gausswin(filtSize)./sum(gausswin(filtSize));

%%
nTrials  = length(clCell);
for n = 1:nTrials
    thistrial = clCell{n};
    [N,edges] = histcounts( thistrial, time );
    Smat(:,n) = nanconv( N', filter);
    Smm(:,n)  = Smat(:,n)./(nTrials*bin);
end;

spikeRate = nanmean(Smat./(nTrials*bin),2);
spikeMat  = Smm';
timeOut   =  linspace( range(1), range(2), length(spikeRate) );

indBase  = find( time >= range(1) & time <= 0 );

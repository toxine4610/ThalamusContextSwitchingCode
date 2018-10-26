
function [spikeRate, spikeMat, time] = makeSpikeRates(foo, range , bin, filtWidth);

clCell = foo;

time = range(1):bin:range(2);
filtSize = round(filtWidth./bin);
filter = gausswin(filtSize)./sum(gausswin(filtSize));

for n = 1:length(clCell)
    test = histc(clCell{n},time);
    if length(test) > 0;
        spikeMat(:,n) = nanconv(histc(clCell{n},time),filter);
    else
        spikeMat(:,n) = 0*time';
    end
end

% spikemMat = nanconv(spikeMat,filter);
%  spikemMat = nanconv(spikeMat,filter');

spikeRate = nanmean(spikeMat')./(bin);

spikeRate = nanconv(spikeRate,filter./norm(filter,1));
spikeMat  = spikeMat';



% error = nanstd(spikeMat)./(length(clCell).*bin);
% 
% errorUp = spikeRate+1.9*mean(error);
% errorDown = spikeRate-1.9*mean(error);
% 
% %errorUP = nanconv(errorUp,filter./sum(filter));
% %errorDown = nanconv(errorDown,filter./sum(filter));
% 
% errorBounds = [errorUp; errorDown];

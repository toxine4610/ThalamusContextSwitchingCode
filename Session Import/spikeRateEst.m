function [spikeRate,error,errorBounds,spikeMat,time] = spikeRateEst(clCell,bin,filtWidth,pre,post)
%clCell= output from mclust? clustered cell? so assumably an array of time
%v. voltage
%bin= length of time bins you're looking at, probably @ sample rate?
%filtWidth=
%pre= ?
%post= last time point

time = -pre:bin:post; %has to be length of clCell
filtSize = round(filtWidth./bin)
% round(X) rounds each element of X to the nearest integer.
    %num bins (samples) per filtwidth, or window for filtering
filter = gausswin(filtSize);
%gausswin(N) returns an N-point Gaussian window.

for n = 1:length(clCell)
    test = histc(clCell{n},time);
    %
    if length(test) > 0;
        spikeMat(:,n) = histc(clCell{n},time);
    else
        spikeMat(:,n) = 0*time';
    end
end

spikemMat = nanconv(spikeMat,filter);
%spikemMat = nanconv(spikeMat,filter');

spikeRate = nanmean(spikeMat')./(bin);
spikeRate = nanconv(spikeRate,filter./sum(filter));
error = nanstd(spikeMat')./(length(clCell).*bin);

% errorUp = spikeRate+1.9*mean(error);
% errorDown = spikeRate-1.9*mean(error);

errorUp = spikeRate+error;
errorDown = spikeRate-error;

%errorUP = nanconv(errorUp,filter./sum(filter));
%errorDown = nanconv(errorDown,filter./sum(filter));

errorBounds = [errorUp errorDown];


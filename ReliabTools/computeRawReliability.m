function [weightedCC, weightedCC_scaled, plottingBins] = computeRawReliability(SpikeTimesCell, maxTime)

% This function plots event synch vs time, and returns scaled for xcorr
% calculation.
% version 1.0, RVR, NYU 2017
% =========================================================================


ct  = 0;
windowType = 'Sliding';
if nargin < 2
maxTime  = 0.6; % 0.6
end;


switch windowType
    case 'Sliding'
        windowLen = 0.05; %0.03
        shiftSize = 0.01;
        binEdges  = 0:shiftSize:(maxTime-windowLen) ;
        plottingBins = binEdges + (windowLen/2);
        nBins = length(binEdges);
    case 'Fixed'
        nBins = 14;
        timeBins = linspace(0,1.5,nBins+1);
end;


P = nchoosek( 1:size(SpikeTimesCell,2), 2);

%%
[ES,ED] = computeReliabTime( nBins, binEdges, windowLen, windowType, SpikeTimesCell, P);

parfor b = 1:nBins
     weightedCC(b) = nanmean(ES{b}).*(size(ES{b},2)./size(P,1)) ;
end;

%%
weightedCC_scaled = Scale_removeNAN(weightedCC);

%%
% [weightedCC_clean, NumPeaks,peakValue, timePeak] = cleanReliabilityTimeSeries( FF_SC, bins_for_plotting, weightedCC, ES, binEdges, windowLen, SpikeTimesCell );

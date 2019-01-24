function [weightedCC_clean, NumPeaks,peakValue, timePeak, plottingBins] = computeReliabilityNew(SpikeTimesCell, FF_SC, bins_for_plotting)

ct  = 0;
windowType = 'Sliding';
maxTime  = 1.5;


switch windowType
    case 'Sliding'
        windowLen = 0.05;
        shiftSize = 0.01;
        binEdges  = 0:shiftSize:(maxTime-windowLen) ;
        plottingBins = binEdges + (windowLen/2);
        nBins = length(binEdges);
    case 'Fixed'
        nBins = 14;
        timeBins = linspace(0,1.5,nBins+1);
end;


%%

P = nchoosek( 1:size(SpikeTimesCell,2), 2);

%%
[ES,ED] = computeReliabTime( nBins, binEdges, windowLen, windowType, SpikeTimesCell, P);

for b = 1:nBins
     weightedCC(b) = nanmean(ES{b}.*size(ES{b},2)./size(P,1)) ;
end;


%%
[weightedCC_clean, NumPeaks,peakValue, timePeak] = cleanReliabilityTimeSeries( FF_SC, bins_for_plotting, weightedCC, ES, binEdges, windowLen, SpikeTimesCell );

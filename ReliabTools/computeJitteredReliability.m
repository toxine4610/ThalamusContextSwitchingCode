

ct  = 0;
windowType = 'Sliding';
maxTime  = 1.5; % 0.6
P = nchoosek( 1:size(SpikeTimesCell,2), 2);


switch windowType
    case 'Sliding'
        windowLen = 0.05; %0.03
        shiftSize = 0.01;
        binEdges  = -2.4:shiftSize:(maxTime-windowLen) ;
        plottingBins = binEdges + (windowLen/2);
        nBins = length(binEdges);
    case 'Fixed'
        nBins = 14;
        timeBins = linspace(0,1.5,nBins+1);
end;


for b = 1 : nBins
    
    switch windowType
        case 'Sliding'
            bin1  = binEdges(b);
            bin2  = binEdges(b)+windowLen;
        case 'Fixed'
            bin1  = timeBins(b);
            bin2  = timeBins(b+1);
    end;
    
    for t = 1:size(SpikeTimesCell,2)
        binnedSpikeTimesCell{b}{t} = SpikeTimesCell{t}( SpikeTimesCell{t} >= bin1 &  SpikeTimesCell{t} <= bin2  );
    end;
    
end;

%%
baselineBins = find( plottingBins >= -0.6 & plottingBins < 0 );
delayBins    = find( plottingBins >= 0    & plottingBins <= 0.6 );

%%
Bins = reshape(binEdges(1:360),[length(delayBins),6]);
for x = 1:size(Bins,2)
    binEdgeslocal = Bins(:,x);
    nBinslocal = length(binEdgeslocal);
    [ES,ED] = computeReliabTime( nBinslocal, binEdgeslocal, windowLen, windowType, SpikeTimesCell, P);
    
    for b = 1:nBinslocal
     weightedCC(b) = nanmean(ES{b}).*(size(ES{b},2)./size(P,1)) ;
    end;

    weightedCC_scaled = Scale_removeNAN(weightedCC);
    
    ShuffledRel{x} = weightedCC;
    clear weightedCC
end;


ci = bootci(500, @(x)nanmean(x,1), y(1:4,:) );
p = anova1( [y(5,:); ci(2,:)]',[], 'off');


%%
% 
% NumRandIter = 20;
% 
% parfor iter = 1:NumRandIter
%     disp(iter);
%     
%     sTrain = [];
%     indBase   = randi( [min(baselineBins), max(baselineBins)], 1 );
%     indDelay  = randi( [min(delayBins), max(delayBins)], 1);
%     
%     indBase = [indBase: indBase+10];
%     indDelay = [indDelay: indDelay+10];
%      
%      
%      
%     base = binnedSpikeTimesCell(indBase);
%     for i = 1:size(base,2)
%         for t = 1:size(SpikeTimesCell,2)
%             base{i}{t} = abs( base{i}{t} );
%         end
%     end
%     
%     sTrain  = binnedSpikeTimesCell;
%     sTrain(indDelay) = base;
%     
% 
%     
%     [weightedCC, weightedCC_scaled] = computeRawReliability_prebinnedSpikeTrains(sTrain, delayBins, nBins, binEdges, windowLen, windowType, P);
%     
%     RelNoise{iter} = weightedCC(delayBins);
%     
%     
% end;

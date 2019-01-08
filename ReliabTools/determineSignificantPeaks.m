function p = determineSignificantPeaks(SpikeTimesCell);

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
baselineBins = find( binEdges >= -1.5 & binEdges < 0 );
delayBins    = find( binEdges >= 0    & binEdges <= 0.5 );

%% OLD method
% Bins = reshape(binEdges(1:360),[length(delayBins),6]);
%
% for x = 5
binEdgeslocal = binEdges(delayBins);
nBinslocal = length(binEdgeslocal);
[ES,ED] = computeReliabTime( nBinslocal, binEdgeslocal, windowLen, windowType, SpikeTimesCell, P);

parfor b = 1:nBinslocal
    weightedCC{b} = nanmean(ES{b}).*(size(ES{b},2)./size(P,1)) ;
end;

weightedCC_scaled = Scale_removeNAN(cell2mat(weightedCC));

%     ShuffledRel{x} = weightedCC;
%     clear weightedCC
% end;

% y = cell2mat(ShuffledRel');
% try
% ci = bootci(500, @(x)nanmean(x,1), y(1:4,:) );
% p = anova1( [y(5,:); ci(2,:)]',[], 'off');
% catch err
%     p  =NaN;
% end
% figure(1);
% plot(ci'); hold on;
% plot(y(5,:),'color','k','linewidth',5);

%% NEW method
% ppm = ParforProgMon('Calculating Null Distribution ', 500 );

parfor iter  = 1:20
    
%     ppm.increment();
    randindex =  randi([1,length(binEdges)-length(delayBins)]);
    binEdgeslocal = binEdges( randindex : randindex+(length(delayBins)-1) );
    nBinslocal = length(binEdgeslocal);
    [ES,ED] = computeReliabTime( nBinslocal, binEdgeslocal, windowLen, windowType, SpikeTimesCell, P);
    
    for b = 1:nBinslocal
        ShuffledRel{iter}(b) = nanmean(ES{b}).*(size(ES{b},2)./size(P,1)) ;
    end;
    %     ShuffledRel{iter} = weightedCC;
end
y = cell2mat(ShuffledRel');

try
    ci = bootci(500, @(x)nanmean(x,1), y );
    p = ranksum( [cell2mat(weightedCC)],  [ci(2,:)],'tail','right' );
catch err
    p  = NaN;
end;


function [param] = computeShuffledReliability(SpikeTimesCell_1, SpikeTimesCell_2);

%%
NumIter = 20;
ModulationIndex = @(x,y) (x-y)./(x+y);

fprintf('Shuffling....\n');

parfor iter = 1:NumIter
%     disp(num2str(iter))
    %% Create 2 shuffled spike trains..
    [Train_1, Test_1] = crossvalind('Resubstitution', length(SpikeTimesCell_1), [0.5,0.5]);
    [Train_2, Test_2] = crossvalind('Resubstitution', length(SpikeTimesCell_2), [0.5,0.5]);
    
%     clear SpikeTimesCell_1T SpikeTimesCell_2T S
    
    SpikeTimesCell_1T = SpikeTimesCell_1(Train_1);
    SpikeTimesCell_2T = SpikeTimesCell_2(Train_2);
    
    S = combineSpikeTrains(SpikeTimesCell_1T, SpikeTimesCell_2T);
    
    S1 = S(randperm(max(size(S))));
    
%     clear SpikeTimesCell_1T SpikeTimesCell_2T S weightedCC_clean
    
    SpikeTimesCell_1T = SpikeTimesCell_1(Test_1);
    SpikeTimesCell_2T = SpikeTimesCell_2(Test_2);
    
    S = combineSpikeTrains(SpikeTimesCell_1T, SpikeTimesCell_2T);
    
    S2 = S(randperm(max(size(S))));
    %% Compute Reliability and similarity between each
    
    [S1_weightedCC, S1_weightedCC_scaled] = computeRawReliability(S1);
    [S2_weightedCC, S2_weightedCC_scaled] = computeRawReliability(S2);
    
    [maxCC{iter}, maxLag{iter}, Consistent{iter}] = computeSimilarityPeaks( S1_weightedCC_scaled, S2_weightedCC_scaled );

    S1_RelTS{iter} = S1_weightedCC;
    S2_RelTS{iter} = S2_weightedCC;
    
end;

EstLag = nanmean( jackknife( @(x) nanmean(x), cell2mat( maxLag)) );
EstCC  = nanmean( jackknife( @(x) nanmean(x), cell2mat( maxCC)) );

param.EstCC         = EstCC;
param.EstLag        = EstLag;
param.maxCC         = maxCC;
param.maxLag        = maxLag;
param.Consistent    = Consistent;
param.S1_RelTS      = S1_RelTS;
param.S2_RelTS      = S2_RelTS;

fprintf('done\n');
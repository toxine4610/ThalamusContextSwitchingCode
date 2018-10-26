

function [PFC, MD, goodPFC, goodMD] = cleanData(Spfc, Smd, Z_C1)

addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/Rajeev_Code/');

%% compile all spikes.....................................................

if isfield( Spfc, 'SpikeTimes_R1C1_inc') == 1
    ver = 'new';
elseif isfield( Spfc, 'SpikeTimes_R1C1_inc') == 0
    ver = 'old';
end;

switch ver
    case 'new'
        % incorrect as well as corrects were saved.
        for i= 1:numel(Spfc)
            
            SpikeTimes = cell(size(Z_C1,1),1);
            
            trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 1;
            
            SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R1C1;
            
            trialSubset = Z_C1(:,1)==0 & Z_C1(:,2) == 1;
            
            SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R1C1_inc;
            
            trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 2;
            
            SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R2C1;
            
            trialSubset = Z_C1(:,1)==0 & Z_C1(:,2) == 2;
            
            SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R2C1_inc;
            
            PFC(i).SpikeTimes = SpikeTimes;
            [spikeRate, spikeMat, timeOut, indBase] = makeSpikeRates(PFC(i).SpikeTimes, [-0.1, 0.7] , 0.001, 0.03);
            PFC(i).muRate = spikeRate;
            PFC(i).Rate   = spikeMat;
            PFC(i).TRRate = nanmean( spikeMat(:, find(timeOut >= 0 & timeOut <= 0.6 )),2);
            PFC(i).Ratenorm = ( spikeRate - nanmin( spikeRate ) ) ./ ( nanmax( spikeRate ) - nanmin( spikeRate ) );
            PFC(i).TRRatenorm = ( PFC(i).TRRate - nanmin(PFC(i).TRRate ) ) ./ ( nanmax( PFC(i).TRRate ) - nanmin( PFC(i).TRRate ) );
            
            [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(PFC(i).SpikeTimes );
            R(i) = findFractionBlankTrials(SpikeCount, bins_for_plotting);
            
        end;
        
        goodPFC = R < 0.55;
        PFC = PFC(goodPFC);
        
    case 'old'
        % no incorrects were saved
        for i= 1:numel(Spfc)
            
            SpikeTimes = cell(size(Z_C1,1),1);
            
            trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 1;
            
            SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R1C1;
            
            trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 2;
            
            SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R2C1;
            
            
            PFC(i).SpikeTimes = SpikeTimes;
            [spikeRate, spikeMat, timeOut, indBase] = makeSpikeRates(PFC(i).SpikeTimes, [-0.1, 0.7] , 0.001, 0.03);
            PFC(i).muRate = spikeRate;
            PFC(i).Rate   = spikeMat;
            PFC(i).TRRate = nanmean( spikeMat(:, find(timeOut >= 0 & timeOut <= 0.6 )),2);
            PFC(i).Ratenorm = ( spikeRate - nanmin( spikeRate ) ) ./ ( nanmax( spikeRate ) - nanmin( spikeRate ) );
            PFC(i).TRRatenorm = ( PFC(i).TRRate - nanmin(PFC(i).TRRate ) ) ./ ( nanmax( PFC(i).TRRate ) - nanmin( PFC(i).TRRate ) );
            
            [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(PFC(i).SpikeTimes );
            R(i) = findFractionBlankTrials(SpikeCount, bins_for_plotting);
            
        end;
        
        goodPFC = R < 0.55;
        PFC = PFC(goodPFC);
end



%%

clear R;

switch ver
    case 'new'
        for i= 1:numel(Smd)
            
            SpikeTimes = cell(size(Z_C1,1),1);
            
            trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 1 ;
            
            SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R1C1;
            
            trialSubset = Z_C1(:,1)==0 & Z_C1(:,2) == 1;
            
            SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R1C1_inc;
            
            trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 2;
            
            SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R2C1;
            
            trialSubset = Z_C1(:,1)==0 & Z_C1(:,2) == 2;
            
            SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R2C1_inc;
            
            MD(i).SpikeTimes = SpikeTimes;
            [spikeRate, spikeMat, timeOut, indBase] = makeSpikeRates(MD(i).SpikeTimes, [-0.1, 0.7] , 0.001, 0.03);
            MD(i).muRate = spikeRate;
            MD(i).Rate   = spikeMat;
            MD(i).TRRate = nanmean( spikeMat(:, find(timeOut >= 0 & timeOut <= 0.6 )),2);
            MD(i).Ratenorm = ( spikeRate - nanmin( spikeRate ) ) ./ ( nanmax( spikeRate ) - nanmin( spikeRate ) );
            MD(i).TRRatenorm = ( MD(i).TRRate - nanmin(MD(i).TRRate ) ) ./ ( nanmax( MD(i).TRRate ) - nanmin(MD(i).TRRate ) );
            
            [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(MD(i).SpikeTimes );
            R(i) = findFractionBlankTrials(SpikeCount, bins_for_plotting);
        end;
        
        goodMD = R < 0.55;
        MD = MD(goodMD);
        
    case 'old'
        for i= 1:numel(Smd)
            
            SpikeTimes = cell(size(Z_C1,1),1);
            
            trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 1 ;
            
            SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R1C1;
            
            trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 2;
            
            SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R2C1;

            MD(i).SpikeTimes = SpikeTimes;
            [spikeRate, spikeMat, timeOut, indBase] = makeSpikeRates(MD(i).SpikeTimes, [-0.1, 0.7] , 0.001, 0.03);
            MD(i).muRate = spikeRate;
            MD(i).Rate   = spikeMat;
            MD(i).TRRate = nanmean( spikeMat(:, find(timeOut >= 0 & timeOut <= 0.6 )),2);
            MD(i).Ratenorm = ( spikeRate - nanmin( spikeRate ) ) ./ ( nanmax( spikeRate ) - nanmin( spikeRate ) );
            MD(i).TRRatenorm = ( MD(i).TRRate - nanmin(MD(i).TRRate ) ) ./ ( nanmax( MD(i).TRRate ) - nanmin(MD(i).TRRate ) );
            
            [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(MD(i).SpikeTimes );
            R(i) = findFractionBlankTrials(SpikeCount, bins_for_plotting);
        end;
        
        goodMD = R < 0.55;
        MD = MD(goodMD);
        
end;
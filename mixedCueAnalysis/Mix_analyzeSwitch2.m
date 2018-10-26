


context = Z_C1(:,9);
dC      = diff(context);
SwitchTrial = find(dC~=0); % switch occurs after this trial

if length(SwitchTrial) == 1
    RangeC1 = 1:SwitchTrial;
    RangeC2 = SwitchTrial+1:size(Z_C1,1);
elseif length(SwitchTrial) == 2
    RangeC1 = 1:SwitchTrial(1);
    RangeC2 = SwitchTrial(1)+1 : SwitchTrial(2);
    RangeC1rep = SwitchTrial(2)+1: size(Z_C1,1);
end;


Z_context1 = Z_C1(RangeC1,:);
Z_context2 = Z_C1(RangeC2,:);
Z_context1rep = Z_C1(RangeC1rep,:);

window1 =  [ SwitchTrial(1) - 10: SwitchTrial(1) + 10];
window2 =  [ SwitchTrial(2) - 10: SwitchTrial(2) + 10];




%% compile all spikes.....................................................

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



clear R;
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


%%

P = nchoosek(1:numel(PFC),2);

for i = 1:length(window1)
    for j = 1:size(P,1)
        clear cell1 cell2
        cell1 = PFC(P(j,1)).Rate(window1(i),find(timeOut >= 0 & timeOut <= 0.6 ));
        cell2 = PFC(P(j,2)).Rate(window1(i),find(timeOut >= 0 & timeOut <= 0.6 ));
        
        CC(i,j) = corr( (cell1)',(cell2)');
    end;
end;
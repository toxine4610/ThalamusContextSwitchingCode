
addpath ../Rajeev_Code
Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/SOMCre/DualModalityCue/';

%%
CMI   = @(X,Y) (X-Y)./(X+Y);
goodT =  @(X) find(~cellfun(@isempty,X));
bin = 0.01;
filtWidth = 0.05;
pre =  0.2;
post = 1.8;

%%
ExptDates = {'2017-12-15','2017-12-16','2017-12-17','2017-12-18','2017-12-20','2017-12-22','2017-12-23','2017-12-24','2017-12-26'};

%%

for d = 1:length(ExptDates)
    
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} '/Somcre_mixed.mat'] )
    
    %%
    
    [Spfc_mix, Smd_mix, ZVis, ZAud] = packageData(Z_C1, Smd, Spfc);
    
    meanPerf_Vis{d} = mean(ZVis(:,8));
    meanPerf_Aud{d} = mean(ZAud(:,8));
    PerfDiff{d}     = (meanPerf_Aud{d} - meanPerf_Vis{d});
    PerfAud{d}      = meanPerf_Aud{d};
    PerfVis{d}      = meanPerf_Vis{d};

    %%
    ct = 0;
    
    for i = 1:numel(Spfc_mix)
        
        fprintf('Processing Cell %d of %d ....', i, numel(Spfc_mix) );
        
        clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateR1C1 spikeRateR2C1 spikeRateR1C2 spikeRateR2C2
        % correct =============================================================
        R1C1  = Spfc_mix(i).SpikeTimes_R1_sound;
        R2C1  = Spfc_mix(i).SpikeTimes_R2_sound;
        
        R1C2  = Spfc_mix(i).SpikeTimes_R1_nosound;
        R2C2  = Spfc_mix(i).SpikeTimes_R2_nosound;
        
        R1C1  =  R1C1( goodT(R1C1) );
        R1C2  =  R1C2( goodT(R1C2) );
        R2C1  =  R2C1( goodT(R2C1) );
        R2C2  =  R2C2( goodT(R2C2) );
        
        
        C1 = [R1C1, R2C1];
        C2 = [R1C2, R2C2];
        
        [mean_spikeRate_1,errorBounds,spikeMat,time] = estimateSpikeRates( C1, 0.01,  0.05, 1, 1.5 );
        [mean_spikeRate_2,errorBounds,spikeMat,time] = estimateSpikeRates( C2, 0.01,  0.05, 1, 1.5 );
        
        [~, CMI_shuffled_delay, ~] = computeShuffledModulation(C1, C2);
        
        base = find( time < 0 );
        delay = find( time >= 0 & time <= 0.2);
        post  = find( time > 0.2 & time <= 1.0 );
        
        CMI_delay(i) = CMI( nanmean( mean_spikeRate_1(delay) ), nanmean( mean_spikeRate_2(delay) ) );
        
        clear perm
        perm = cell2mat( CMI_shuffled_delay );
        p(i) = sum(abs(perm) > abs(CMI_delay(i)))/100;
        
        fprintf('... Done\n');
    end;
    CMI_PCD{d} = CMI_delay;
    Prob{d}    = p;
    clear CMI_delay p
end;

%%
for d = 1:length(ExptDates)
    
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} '/Somcre_mixed.mat'] )
    
    %%
    
    [Spfc_mix, Smd_mix, ZVis, ZAud] = packageData(Z_C1, Smd, Spfc);
    
    meanPerf_Vis{d} = mean(ZVis(:,8));
    meanPerf_Aud{d} = mean(ZAud(:,8));
    PerfDiff{d}     = (meanPerf_Aud{d} - meanPerf_Vis{d});
    PerfAud{d}      = meanPerf_Aud{d};
    PerfVis{d}      = meanPerf_Vis{d};

    %%
    ct = 0;

    for i = 1:numel(Smd_mix)
        
        fprintf('Processing Cell %d of %d ....', i, numel(Smd_mix) );
        
        clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateR1C1 spikeRateR2C1 spikeRateR1C2 spikeRateR2C2
        % correct =============================================================
        R1C1  = Smd_mix(i).SpikeTimes_R1_sound;
        R2C1  = Smd_mix(i).SpikeTimes_R2_sound;
        
        R1C2  = Smd_mix(i).SpikeTimes_R1_nosound;
        R2C2  = Smd_mix(i).SpikeTimes_R2_nosound;
        
        R1C1  =  R1C1( goodT(R1C1) );
        R1C2  =  R1C2( goodT(R1C2) );
        R2C1  =  R2C1( goodT(R2C1) );
        R2C2  =  R2C2( goodT(R2C2) );
        
        
        C1 = [R1C1, R2C1];
        C2 = [R1C2, R2C2];
        
        [mean_spikeRate_1,errorBounds,spikeMat,time] = estimateSpikeRates( C1, 0.01,  0.05, 1, 1.5 );
        [mean_spikeRate_2,errorBounds,spikeMat,time] = estimateSpikeRates( C2, 0.01,  0.05, 1, 1.5 );
        
        [~, CMI_shuffled_delay, ~] = computeShuffledModulation(C1, C2);
        
        base = find( time < 0 );
        delay = find( time >= 0 & time <= 0.2);
        post  = find( time > 0.2 & time <= 1.0 );
        
        MD_CMI_delay(i) = CMI( nanmean( mean_spikeRate_1(delay) ), nanmean( mean_spikeRate_2(delay) ) );
        
        clear perm
        perm = cell2mat( CMI_shuffled_delay );
        MD_p(i) = sum(abs(perm) > abs(MD_CMI_delay(i)))/100;
        
        fprintf('... Done\n');
    end;
    MD_CMI_PCD{d} = MD_CMI_delay;
    MD_Prob{d}    = MD_p;
    clear CMI_delay p
    
end;

%%
cmiPcd = cell2mat(MD_CMI_PCD);
pval   = cell2mat(MD_Prob);
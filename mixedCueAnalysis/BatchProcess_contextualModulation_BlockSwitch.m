

function BatchProcess_contextualModulation_BlockSwitch


addpath('../Rajeev_Code/')


%% Compile all the cells....
warning off;
datapath = '../ForContextSwitchProject/DoubleCueDatabase/Production Learning MDterminal';

X = dir(datapath);
MouseID = {X(3:end).name};
MyMouse  = { 'BPM13', 'BPM25', 'BPM43' };

SmdMega_av = [];
SpfcMega_av = [];

SmdMega_va = [];
SpfcMega_va = [];


ct = 0;

for m = 1:3
    mouse = MouseID{m};
    fprintf('Processing Mouse: %s\n', mouse );
    dates_this_mouse = getDateList([datapath filesep mouse]);
    for d = 1:length(dates_this_mouse)
        ct = ct+1;
        fprintf('  Processing date: %s\n', dates_this_mouse{d});
        ST = load( [datapath filesep mouse filesep dates_this_mouse{d} filesep MyMouse{m} '_learning_mdterminal_session' ] );
        [Spfc_mix, Smd_mix, ZVis, ZAud, first] = packageData_BlockSwitch_laser(ST.Z_C1, ST.Smd, ST.Spfc);
        [ExpH, perf1, Toss] = getBehavior(ST.Z_C1);
        
        returnBehaviorPlot( ST.Z_C1, ct)
        
        Context(ct)  = first;
        avePerf(ct) = median(ExpH);
        pCorrect(ct) = length( find(Toss == 1) )./length( Toss );
        
        if Context(ct) == 0
            SpfcMega_av = [SpfcMega_av, Spfc_mix];
            SmdMega_av  = [SmdMega_av, Smd_mix];
        elseif Context(ct) == 3
            SpfcMega_va = [SpfcMega_va, Spfc_mix];
            SmdMega_va  = [SmdMega_va, Smd_mix];
        end;
        
        clear Spfc_mix Smd_mix
    end;
end;

%% Spike timing parameters .........

CMI   = @(X,Y) (X-Y)./(X+Y);
goodT =  @(X) find(~cellfun(@isempty,X));
bin = 0.01;
filtWidth = 0.05;
pre =  0.2;
post = 1.5;


%% PFC ====================================================================
brokenCells = [];
ct = 0;

for i = 1:numel(SpfcMega_av)
    fprintf('Processing PFC Cell %d\n', i );
    clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateC1 spikeRateC2 spikeMatC1 spikeMatC2
    
    % correct, no laser ===================================================
    R1C1  = SpfcMega_av(i).SpikeTimes_R1_sound_NoLaser_first;
    R2C1  = SpfcMega_av(i).SpikeTimes_R2_sound_NoLaser_first;
    
    R1C2  = SpfcMega_av(i).SpikeTimes_R1_nosound_NoLaser;
    R2C2  = SpfcMega_av(i).SpikeTimes_R2_nosound_NoLaser;
    
    R1C3  = SpfcMega_av(i).SpikeTimes_R1_sound_NoLaser_rep;
    R2C3  = SpfcMega_av(i).SpikeTimes_R2_sound_NoLaser_rep;
    
    
    R1C1  =  R1C1( goodT(R1C1) );
    R1C2  =  R1C2( goodT(R1C2) );
    R1C3  =  R1C3( goodT(R1C3) );
    
    R2C1  =  R2C1( goodT(R2C1) );
    R2C2  =  R2C2( goodT(R2C2) );
    R2C3  =  R2C3( goodT(R2C3) );
    
    C1 = [R1C1, R2C1];
    C2 = [R1C2, R2C2];
    C3 = [R1C2, R2C3];
    
    [spikeRateC1,errorBounds,spikeMatC1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
    [spikeRateC2,errorBounds,spikeMatC2,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
    [spikeRateC3,errorBounds,spikeMatC3,time] = estimateSpikeRates( C3, bin, filtWidth, pre, post );
    
    indBase = find( time >= -0.25 & time  <=0 );
    indDelay = find( time >= 0.0 & time <= 0.6 );
    
    PFC_muFR_CFD_C1(i) = nanmean(spikeRateC1(indDelay));
    PFC_muFR_CFD_C2(i) = nanmean(spikeRateC2(indDelay));
    PFC_muFR_CFD_C3(i) = nanmean(spikeRateC3(indDelay));
    
    PFC_muFR_BAS_C1(i) = nanmean(spikeRateC1(indBase));
    PFC_muFR_BAS_C2(i) = nanmean(spikeRateC2(indBase));
    PFC_muFR_BAS_C3(i) = nanmean(spikeRateC3(indBase));
    
    % Context 1-->2
    PFC_CMI_BAS(i)  = CMI( PFC_muFR_BAS_C1(i), PFC_muFR_BAS_C2(i) );
    PFC_CMI_CFD(i)  = CMI( PFC_muFR_CFD_C1(i), PFC_muFR_CFD_C2(i) );
    
    % Context 2-->3
    PFC_CMI_BAS_2nd(i)  = CMI( PFC_muFR_BAS_C2(i), PFC_muFR_BAS_C3(i) );
    PFC_CMI_CFD_2nd(i)  = CMI( PFC_muFR_CFD_C2(i), PFC_muFR_CFD_C3(i) );
    
    % Context 1-->3
    PFC_CMI_BAS_3nd(i)  = CMI( PFC_muFR_BAS_C1(i), PFC_muFR_BAS_C3(i) );
    PFC_CMI_CFD_3nd(i)  = CMI( PFC_muFR_CFD_C1(i), PFC_muFR_CFD_C3(i) );
    
    PFC_CMI_C1_we(i)  =  CMI( PFC_muFR_CFD_C1(i), PFC_muFR_BAS_C1(i));
    PFC_CMI_C2_we(i)  =  CMI( PFC_muFR_CFD_C2(i), PFC_muFR_BAS_C2(i));
    PFC_CMI_C3_we(i)  =  CMI( PFC_muFR_CFD_C3(i), PFC_muFR_BAS_C3(i));
    
    PFC_context_dp(i,:) = (nanmean( spikeMatC1' ) - nanmean( spikeMatC2' ))./ sqrt( 0.5.*( var( spikeMatC1',[],1) + var( spikeMatC2',[], 1) ) );
    
    
    clear C1 C2 C3
end;


%% MD ====================================================================
brokenCells = [];
ct = 0;

for i = 1:numel(SmdMega_av)
    fprintf('Processing MD Cell %d\n', i );
    clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateC1 spikeRateC2 spikeMatC1 spikeMatC2
    
    % correct, no laser ===================================================
    R1C1  = SmdMega_av(i).SpikeTimes_R1_sound_NoLaser_first;
    R2C1  = SmdMega_av(i).SpikeTimes_R2_sound_NoLaser_first;
    
    R1C2  = SmdMega_av(i).SpikeTimes_R1_nosound_NoLaser;
    R2C2  = SmdMega_av(i).SpikeTimes_R2_nosound_NoLaser;
    
    R1C3  = SmdMega_av(i).SpikeTimes_R1_sound_NoLaser_rep;
    R2C3  = SmdMega_av(i).SpikeTimes_R2_sound_NoLaser_rep;
    
    
    R1C1  =  R1C1( goodT(R1C1) );
    R1C2  =  R1C2( goodT(R1C2) );
    R1C3  =  R1C3( goodT(R1C3) );
    
    R2C1  =  R2C1( goodT(R2C1) );
    R2C2  =  R2C2( goodT(R2C2) );
    R2C3  =  R2C3( goodT(R2C3) );
    
    C1 = [R1C1, R2C1];
    C2 = [R1C2, R2C2];
    C3 = [R1C2, R2C3];
    
    [spikeRateC1,errorBounds,spikeMatC1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
    [spikeRateC2,errorBounds,spikeMatC2,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
    [spikeRateC3,errorBounds,spikeMatC3,time] = estimateSpikeRates( C3, bin, filtWidth, pre, post );
    
    indBase = find( time >= -0.25 & time  <=0 );
    indDelay = find( time >= 0.0 & time <= 0.6 );
    
    MD_muFR_CFD_C1(i) = nanmean(spikeRateC1(indDelay));
    MD_muFR_CFD_C2(i) = nanmean(spikeRateC2(indDelay));
    MD_muFR_CFD_C3(i) = nanmean(spikeRateC3(indDelay));
    
    MD_muFR_BAS_C1(i) = nanmean(spikeRateC1(indBase));
    MD_muFR_BAS_C2(i) = nanmean(spikeRateC2(indBase));
    MD_muFR_BAS_C3(i) = nanmean(spikeRateC3(indBase));
    
    % Context 1-->2
    MD_CMI_BAS(i)  = CMI( MD_muFR_BAS_C1(i), MD_muFR_BAS_C2(i) );
    MD_CMI_CFD(i)  = CMI( MD_muFR_CFD_C1(i), MD_muFR_CFD_C2(i) );
    
    % Context 2-->3
    MD_CMI_BAS_2nd(i)  = CMI( MD_muFR_BAS_C2(i), MD_muFR_BAS_C3(i) );
    MD_CMI_CFD_2nd(i)  = CMI( MD_muFR_CFD_C2(i), MD_muFR_CFD_C3(i) );
    
    % Context 1-->3
    MD_CMI_BAS_3nd(i)  = CMI( MD_muFR_BAS_C1(i), MD_muFR_BAS_C3(i) );
    MD_CMI_CFD_3nd(i)  = CMI( MD_muFR_CFD_C1(i), MD_muFR_CFD_C3(i) );
    
    MD_CMI_C1_we(i)  =  CMI( MD_muFR_CFD_C1(i), MD_muFR_BAS_C1(i));
    MD_CMI_C2_we(i)  =  CMI( MD_muFR_CFD_C2(i), MD_muFR_BAS_C2(i));
    MD_CMI_C3_we(i)  =  CMI( MD_muFR_CFD_C3(i), MD_muFR_BAS_C3(i));
    
    MD_context_dp(i,:) = (nanmean( spikeMatC1' ) - nanmean( spikeMatC2' ))./ sqrt( 0.5.*( var( spikeMatC1',[],1) + var( spikeMatC2',[], 1) ) );

    clear C1 C2 C3
end;


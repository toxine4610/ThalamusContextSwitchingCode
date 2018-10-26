
% function BatchProcess_contextualModulation_withLaser


addpath('../Rajeev_Code/')


%% Compile all the cells....

datapath = '../ForContextSwitchProject/DoubleCueDatabase/Production PFCTerminal Unilat';

X = dir(datapath);
MouseID = {X(3:end).name};
MyMouse  = { 'BPM13', 'BPM25', 'BPM43' };

SmdMega = [];
SpfcMega = [];

ct = 0;

for m = 1:length(MouseID)
    mouse = MouseID{m};
    fprintf('Processing Mouse: %s\n', mouse );
    dates_this_mouse = getDateList([datapath filesep mouse]);
    for d = 1:length(dates_this_mouse)
        ct = ct+1;
        fprintf('  Processing date: %s\n', dates_this_mouse{d});
        ST = load( [datapath filesep mouse filesep dates_this_mouse{d} filesep MyMouse{m} '_pfcterminal_session' ] );
        [Spfc_mix, Smd_mix, ZVis, ZAud] = packageData_Laser(ST.Z_C1, ST.Smd, ST.Spfc);
        [ExpH, perf1, Toss] = getBehavior(ST.Z_C1);
        
        avePerf(ct) = median(ExpH);
        pCorrect(ct) = length( find(Toss == 1) )./length( Toss );
      
        SpfcMega = [SpfcMega, Spfc_mix];
        SmdMega  = [SmdMega, Smd_mix];
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

for i = 1:numel(SpfcMega)
    fprintf('Processing PFC Cell %d\n', i );
    
    try
        clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateC1 spikeRateC2 spikeMatC1 spikeMatC2
        
        % correct, no laser ===================================================
        R1C1  = SpfcMega(i).SpikeTimes_R1_sound_NoLaser;
        R2C1  = SpfcMega(i).SpikeTimes_R2_sound_NoLaser;
        
        R1C2  = SpfcMega(i).SpikeTimes_R1_nosound_NoLaser;
        R2C2  = SpfcMega(i).SpikeTimes_R2_nosound_NoLaser;
        
        R1C1  =  R1C1( goodT(R1C1) );
        R1C2  =  R1C2( goodT(R1C2) );
        R2C1  =  R2C1( goodT(R2C1) );
        R2C2  =  R2C2( goodT(R2C2) );
        
        C1 = [R1C1, R2C1];
        C2 = [R1C2, R2C2];
        
        [spikeRateC1,errorBounds,spikeMatC1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
        [spikeRateC2,errorBounds,spikeMatC2,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
        
        indBase = find( time >= -0.25 & time  <=0 );
        indDelay = find( time >= 0.0 & time <= 0.6 );
        
        PFC_muFR_CFD_C1(i) = nanmean(spikeRateC1(indDelay));
        PFC_muFR_CFD_C2(i) = nanmean(spikeRateC2(indDelay));
        
        PFC_muFR_BAS_C1(i) = nanmean(spikeRateC1(indBase));
        PFC_muFR_BAS_C2(i) = nanmean(spikeRateC2(indBase));
        
        PFC_CMI_BAS(i)  = CMI( PFC_muFR_BAS_C1(i), PFC_muFR_BAS_C2(i));
        PFC_CMI_CFD(i)  = CMI( PFC_muFR_CFD_C1(i), PFC_muFR_CFD_C2(i));
        
        PFC_CMI_C1_we(i)  =  CMI( PFC_muFR_CFD_C1(i), PFC_muFR_BAS_C1(i));
        PFC_CMI_C2_we(i)  =  CMI( PFC_muFR_CFD_C2(i), PFC_muFR_BAS_C2(i));
        
        
        PFC_context_dp(i,:) = (nanmean( spikeMatC1' ) - nanmean( spikeMatC2' ))./ sqrt( 0.5.*( var( spikeMatC1',[],1) + var( spikeMatC2',[], 1) ) );
        
        
        clear C1 C2
        
        C1 = [R1C1, R2C2]; % (HP-->Vis, Green --> Aud)
        C2 = [R1C2, R2C1]; % (UV-->Vis, LP-->Aud)
        
        [spikeRateC1_bl,errorBounds,spikeMatC1_bl,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
        [spikeRateC2_bl,errorBounds,spikeMatC2_bl,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
        
        indBase = find( time >= -0.25 & time  <=0 );
        indDelay = find( time >= 0.0 & time <= 0.6 );
        
        PFC_muFR_CFD_C1_bl(i) = nanmean(spikeRateC1_bl(indDelay));
        PFC_muFR_CFD_C2_bl(i) = nanmean(spikeRateC2_bl(indDelay));
        
        PFC_muFR_BAS_C1_bl(i) = nanmean(spikeRateC1_bl(indBase));
        PFC_muFR_BAS_C2_bl(i) = nanmean(spikeRateC2_bl(indBase));
        
        PFC_CMI_BAS_bl(i)  = CMI( PFC_muFR_BAS_C1_bl(i), PFC_muFR_BAS_C2_bl(i));
        PFC_CMI_CFD_bl(i)  = CMI( PFC_muFR_CFD_C1_bl(i), PFC_muFR_CFD_C2_bl(i));
        
        PFC_CMI_C1_we_bl(i)  =  CMI( PFC_muFR_CFD_C1_bl(i), PFC_muFR_BAS_C1_bl(i));
        PFC_CMI_C2_we_bl(i)  =  CMI( PFC_muFR_CFD_C2_bl(i), PFC_muFR_BAS_C2_bl(i));
        
        
        P_context_bl_dp(i,:) = (nanmean( spikeMatC1_bl' ) - nanmean( spikeMatC2_bl' ))./ sqrt( 0.5.*( var( spikeMatC1_bl',[],1) + var( spikeMatC2_bl',[], 1) ) );
        
        
        
        clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateC1 spikeRateC2 spikeMatC1 spikeMatC2
        % correct, Laser ======================================================
        
        R1C1  = SpfcMega(i).SpikeTimes_R1_sound_Laser;
        R2C1  = SpfcMega(i).SpikeTimes_R2_sound_Laser;
        
        R1C2  = SpfcMega(i).SpikeTimes_R1_nosound_Laser;
        R2C2  = SpfcMega(i).SpikeTimes_R2_nosound_Laser;
        
        R1C1  =  R1C1( goodT(R1C1) );
        R1C2  =  R1C2( goodT(R1C2) );
        R2C1  =  R2C1( goodT(R2C1) );
        R2C2  =  R2C2( goodT(R2C2) );
        
        C1 = [R1C1, R2C1];
        C2 = [R1C2, R2C2];
        
        [spikeRateC1err,errorBounds,spikeMatC1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
        [spikeRateC2err,errorBounds,spikeMatC2,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
        
        indBase = find( time >= -0.25 & time  <=0 );
        indDelay = find( time >= 0.0 & time <= 0.6 );
        
        PFC_muFR_CFD_C1err(i) = nanmean(spikeRateC1err(indDelay));
        PFC_muFR_CFD_C2err(i) = nanmean(spikeRateC2err(indDelay));
        
        PFC_muFR_BAS_C1err(i) = nanmean(spikeRateC1err(indBase));
        PFC_muFR_BAS_C2err(i) = nanmean(spikeRateC2err(indBase));
        
        PFC_CMI_BASerr(i)  = CMI( PFC_muFR_BAS_C1err(i), PFC_muFR_BAS_C2err(i));
        PFC_CMI_CFDerr(i)  = CMI( PFC_muFR_CFD_C1err(i), PFC_muFR_CFD_C2err(i));
        
        PFC_CMI_C1_weerr(i)  =  CMI( PFC_muFR_CFD_C1err(i), PFC_muFR_BAS_C1err(i));
        PFC_CMI_C2_weerr(i)  =  CMI( PFC_muFR_CFD_C2err(i), PFC_muFR_BAS_C2err(i));
        
        PFC_context_err_dp(i,:) = (nanmean( spikeMatC1' ) - nanmean( spikeMatC2' ))./ sqrt( 0.5.*( var( spikeMatC1',[],1) + var( spikeMatC2',[], 1) ) );
    catch err
        brokenCells = [brokenCells, i];
    end
    
end;

%%
% for i = 1:numel(SpfcMega)
%     foo = smooth ( PFC_context_dp(i, : ) );
%     DP(i,:) = foo;
%     [v(i), ind(i)] = max( foo );
%
%     foo2 = smooth( PFC_context_err_dp(i,:) );
%     DPerr(i,:) = foo2;
% end;
%
% [~, index] = sort( ind );
% figure(1); subplot(1,2,1);
% imagesc( DP( index, :) );
% axis square; box off;
%
% figure(1); subplot(1,2,2);
% imagesc( DPerr( index, :) );
% axis square; box off



%% MD  ====================================================================

brokenMDCells = [] ;
ct = 0;

for i = 1:numel(SmdMega)
    
    clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateC1 spikeRateC2 spikeMatC1 spikeMatC2
    
    try
    % correct =============================================================
    R1C1  = SmdMega(i).SpikeTimes_R1_sound_NoLaser;
    R2C1  = SmdMega(i).SpikeTimes_R2_sound_NoLaser;
    
    R1C2  = SmdMega(i).SpikeTimes_R1_nosound_NoLaser;
    R2C2  = SmdMega(i).SpikeTimes_R2_nosound_NoLaser;
    
    R1C1  =  R1C1( goodT(R1C1) );
    R1C2  =  R1C2( goodT(R1C2) );
    R2C1  =  R2C1( goodT(R2C1) );
    R2C2  =  R2C2( goodT(R2C2) );
    
    C1 = [R1C1, R2C1];
    C2 = [R1C2, R2C2];
    
    [spikeRateC1,errorBounds,spikeMatC1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
    [spikeRateC2,errorBounds,spikeMatC2,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
    
    indBase = find( time >= -0.25 & time  <=0 );
    indDelay = find( time >= 0.0 & time <= 0.6 );
    
    MD_muFR_CFD_C1(i) = nanmean(spikeRateC1(indDelay));
    MD_muFR_CFD_C2(i) = nanmean(spikeRateC2(indDelay));
    
    MD_muFR_BAS_C1(i) = nanmean(spikeRateC1(indBase));
    MD_muFR_BAS_C2(i) = nanmean(spikeRateC2(indBase));
    
    MD_CMI_BAS(i)  = CMI( MD_muFR_BAS_C1(i), MD_muFR_BAS_C2(i));
    MD_CMI_CFD(i)  = CMI( MD_muFR_CFD_C1(i), MD_muFR_CFD_C2(i));
    
    MD_CMI_C1_we(i)  =  CMI( MD_muFR_CFD_C1(i), MD_muFR_BAS_C1(i));
    MD_CMI_C2_we(i)  =  CMI( MD_muFR_CFD_C2(i), MD_muFR_BAS_C2(i));
    
    
    
    clear C1 C2
    
    C1 = [R1C1, R2C2];
    C2 = [R1C2, R2C1];
    
    [spikeRateC1_bl,errorBounds,spikeMatC1_bl,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
    [spikeRateC2_bl,errorBounds,spikeMatC2_bl,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
    
    indBase = find( time >= -0.25 & time  <=0 );
    indDelay = find( time >= 0.0 & time <= 0.6 );
    
    MD_muFR_CFD_C1_bl(i) = nanmean(spikeRateC1_bl(indDelay));
    MD_muFR_CFD_C2_bl(i) = nanmean(spikeRateC2_bl(indDelay));
    
    MD_muFR_BAS_C1_bl(i) = nanmean(spikeRateC1_bl(indBase));
    MD_muFR_BAS_C2_bl(i) = nanmean(spikeRateC2_bl(indBase));
    
    MD_CMI_BAS_bl(i)  = CMI( MD_muFR_BAS_C1_bl(i), MD_muFR_BAS_C2_bl(i));
    MD_CMI_CFD_bl(i)  = CMI( MD_muFR_CFD_C1_bl(i), MD_muFR_CFD_C2_bl(i));
    
    MD_CMI_C1_we_bl(i)  =  CMI( MD_muFR_CFD_C1_bl(i), MD_muFR_BAS_C1_bl(i));
    MD_CMI_C2_we_bl(i)  =  CMI( MD_muFR_CFD_C2_bl(i), MD_muFR_BAS_C2_bl(i));
    
    
    MD_context_bl_dp(i,:) = (nanmean( spikeMatC1_bl' ) - nanmean( spikeMatC2_bl' ))./ sqrt( 0.5.*( var( spikeMatC1_bl',[],1) + var( spikeMatC2_bl',[], 1) ) );
    
    
    clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateC1 spikeRateC2 spikeMatC1 spikeMatC2
    % errors ============================================================
    
    R1C1  = SmdMega(i).SpikeTimes_R1_sound_Laser;
    R2C1  = SmdMega(i).SpikeTimes_R2_sound_Laser;
    
    R1C2  = SmdMega(i).SpikeTimes_R1_nosound_Laser;
    R2C2  = SmdMega(i).SpikeTimes_R2_nosound_Laser;
    
    R1C1  =  R1C1( goodT(R1C1) );
    R1C2  =  R1C2( goodT(R1C2) );
    R2C1  =  R2C1( goodT(R2C1) );
    R2C2  =  R2C2( goodT(R2C2) );
    
    C1 = [R1C1, R2C1];
    C2 = [R1C2, R2C2];
    
    [spikeRateC1err,errorBounds,spikeMatC1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
    [spikeRateC2err,errorBounds,spikeMatC2,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
    
    indBase = find( time >= -0.25 & time  <=0 );
    indDelay = find( time >= 0.0 & time <= 0.6 );
    
    MD_muFR_CFD_C1err(i) = nanmean(spikeRateC1err(indDelay));
    MD_muFR_CFD_C2err(i) = nanmean(spikeRateC2err(indDelay));
    
    MD_muFR_BAS_C1err(i) = nanmean(spikeRateC1err(indBase));
    MD_muFR_BAS_C2err(i) = nanmean(spikeRateC2err(indBase));
    
    MD_CMI_BASerr(i)  = CMI( MD_muFR_BAS_C1err(i), MD_muFR_BAS_C2err(i));
    MD_CMI_CFDerr(i)  = CMI( MD_muFR_CFD_C1err(i), MD_muFR_CFD_C2err(i));
    
    MD_CMI_C1_weerr(i)  =  CMI( MD_muFR_CFD_C1err(i), MD_muFR_BAS_C1err(i));
    MD_CMI_C2_weerr(i)  =  CMI( MD_muFR_CFD_C2err(i), MD_muFR_BAS_C2err(i));
    
    MD_context_err_dp(i,:) = (nanmean( spikeMatC1' ) - nanmean( spikeMatC2' ))./ sqrt( 0.5.*( var( spikeMatC1',[],1) + var( spikeMatC2',[], 1) ) );
    catch err
        brokenMDCells = [brokenMDCells,i];
    end
end;


%% Make figure

figure(1); set(gcf,'color','w');
% plot( abs(MD_CMI_BAS), abs(MD_CMI_BASerr) ,'o','markersize',10,'markeredgecolor','w','markerfacecolor','g' ); hold on;
plot( abs(PFC_CMI_CFD), abs(PFC_CMI_CFDerr) ,'o','markersize',10,'markeredgecolor','w','markerfacecolor','b' ); hold on;
plot( linspace(0,1,2000), linspace(0,1,2000),'r');

D1 = (abs(MD_CMI_CFDerr) ) - (abs(MD_CMI_CFD));
D2 =  (abs(MD_CMI_BASerr) -  abs(MD_CMI_BAS));

figure(2); set(gcf,'color','w');
histogram( [D1,D2], linspace(-1,1,25), 'normalization','probability'); hold on
histogram( D2, linspace(-1,1,25), 'normalization','probability'); hold on


% figure(2); set(gcf,'color','w');
% plot(  abs(PFC_CMI_BAS), abs(PFC_CMI_BAS_bl) ,'o','markersize',10,'markeredgecolor','w','markerfacecolor','k' ); hold on;
% plot(  abs(PFC_CMI_CFD), abs(PFC_CMI_CFD_bl) ,'o','markersize',10,'markeredgecolor','w','markerfacecolor','k' ); hold on;
% 

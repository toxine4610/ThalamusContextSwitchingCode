
addpath ../Rajeev_Code
% Datapath =  'C:\Users\Halassalab-CG\Dropbox (Personal)\Rajeev\ForContextSwitchProject\DoubleCueDatabase\BPM3_3\CatchTrials\';
% Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/SOMCre/ContextSwitchSSFO/';

%%
CMI   =  @(X,Y) (X-Y)./(X+Y);
goodT =  @(X) find(~cellfun(@isempty,X));
bin = 0.01;
filtWidth = 0.05;
pre =  0.2;
post = 1.8;

%%
%  ExptDates = {'2018-02-11','2018-02-12','2018-02-13', '2018-02-14', '2018-02-16'};
% ExptDates = {'2018-02-17','2018-02-18','2018-02-20','2018-02-22','2018-02-23','2018-02-24'};
% ExptDates = {'2018-03-18', '2018-03-19','2018-03-20','2018-03-21','2018-03-22','2018-03-24'};

% ExptDates = {'2018-04-08', '2018-04-09','2018-04-10'};

Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/BPM3_3/CatchTrials/';
ExptDates = {'2018-04-08', '2018-04-09', '2018-04-10'};
%%
Vmega = [];

for d = 1:length(ExptDates)
    
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} '/BPM33_catchTrials_Session.mat'] )
%     load( [Datapath ExptDates{d} '/Somcre_mixed_ssfo_C1.mat'] )
    Vmega = cat(1, Vmega, Z_C1 );
    %%
    
    [Spfc_mix, Smd_mix, ZVis, ZAud, first, ErrorFrac] = packageData_Laser(Z_C1, Smd, Spfc);
    FirstContext{d} = first;
    EF{d} = ErrorFrac;
    %%
    ct = 0;
    for i = 1:numel(Spfc_mix)
        
        try
            
            clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateR1C1 spikeRateR2C1 spikeRateR1C2 spikeRateR2C2
            % correct nolaser =====================================================
            R1C1  = Spfc_mix(i).SpikeTimes_R1_sound_NoLaser;
            R2C1  = Spfc_mix(i).SpikeTimes_R2_sound_NoLaser;
            
            R1C2  = Spfc_mix(i).SpikeTimes_R1_nosound_NoLaser;
            R2C2  = Spfc_mix(i).SpikeTimes_R2_nosound_NoLaser;
            
            R1C1  =  R1C1( goodT(R1C1) );
            R1C2  =  R1C2( goodT(R1C2) );
            R2C1  =  R2C1( goodT(R2C1) );
            R2C2  =  R2C2( goodT(R2C2) );
            
            
            C1 = [R1C1, R2C1];
            C2 = [R1C2, R2C2];
            
            [spikeRateC1,errorBounds,spikeMatC1Las,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
            [spikeRateC2,errorBounds,spikeMatC2Las,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
            
            % == Dprime within context.
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time2] = estimateSpikeRates( R1C1, bin, filtWidth, 0.1, 0.9 );
            [spikeRateR1C2,errorBounds,spikeMatR1C2,time2] = estimateSpikeRates( R1C2, bin, filtWidth, 0.1, 0.9 );
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time2] = estimateSpikeRates( R2C1, bin, filtWidth, 0.1, 0.9 );
            [spikeRateR2C2,errorBounds,spikeMatR2C2,time2] = estimateSpikeRates( R2C2, bin, filtWidth, 0.1, 0.9 );
            
            DprimeC1{d}(i,:) = ( nanmean( spikeMatR1C1,2 )- nanmean( spikeMatR2C1,2) )./ sqrt( nanvar(spikeMatR1C1,[],2) + nanvar(spikeMatR2C1,[],2) );
            DprimeC2{d}(i,:) = ( nanmean( spikeMatR1C2,2 )- nanmean( spikeMatR2C2,2) )./ sqrt( nanvar(spikeMatR1C2,[],2) + nanvar(spikeMatR2C2,[],2) );

            
            ST = [C1, C2 ];
            [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(ST);
            base = mean_SC( find(bins_for_plotting < 0) );
            resp = mean_SC( find(bins_for_plotting > 0) );
            p{d}(i) = ranksum( base, resp );
            
            FracBlank{d}(i) = findFractionBlankTrials(SpikeCount, bins_for_plotting);
            
            
            % correct laser =====================================================
            R1C1Las  = Spfc_mix(i).SpikeTimes_R1_sound_Laser;
            R2C1Las  = Spfc_mix(i).SpikeTimes_R2_sound_Laser;
            
            R1C2Las  = Spfc_mix(i).SpikeTimes_R1_nosound_Laser;
            R2C2Las  = Spfc_mix(i).SpikeTimes_R2_nosound_Laser;
            
            R1C1Las  =  R1C1Las;
            R1C2Las  =  R1C2Las;
            R2C1Las  =  R2C1Las;
            R2C2Las  =  R2C2Las;
            
            C1Las = [R1C1Las, R2C1Las];
            C2Las = [R1C2Las, R2C2Las];
            
            ST = [C1Las, C2Las ];
            [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(ST);
            base = mean_SC( find(bins_for_plotting < 0) );
            resp = mean_SC( find(bins_for_plotting > 0) );
            pLas{d}(i) = ranksum( base, resp );
            FracBlankLaser{d}(i) = findFractionBlankTrials(SpikeCount, bins_for_plotting);

            
            [spikeRateC1Las,errorBounds,spikeMatC1Las,time] = estimateSpikeRates( C1Las, bin, filtWidth, pre, post );
            [spikeRateC2Las,errorBounds,spikeMatC2Las,time] = estimateSpikeRates( C2Las, bin, filtWidth, pre, post );
            
             % == Dprime within context.
            [spikeRateR1C1,errorBounds,spikeMatR1C1Las,time2] = estimateSpikeRates( R1C1Las, bin, filtWidth, pre, post );
            [spikeRateR1C2,errorBounds,spikeMatR1C2Las,time2] = estimateSpikeRates( R1C2Las, bin, filtWidth, pre, post );
            [spikeRateR2C1,errorBounds,spikeMatR2C1Las,time2] = estimateSpikeRates( R2C1Las, bin, filtWidth, pre, post );
            [spikeRateR2C2,errorBounds,spikeMatR2C2Las,time2] = estimateSpikeRates( R2C2Las, bin, filtWidth, pre, post );
            
            DprimeC1Las{d}(i,:) = ( nanmean( spikeMatR1C1Las,2 )- nanmean( spikeMatR2C1Las,2) )./ sqrt( nanvar(spikeMatR1C1Las,[],2) + nanvar(spikeMatR2C1Las,[],2) );
            DprimeC2Las{d}(i,:) = ( nanmean( spikeMatR1C2Las,2 )- nanmean( spikeMatR2C2Las,2) )./ sqrt( nanvar(spikeMatR1C2Las,[],2) + nanvar(spikeMatR2C2Las,[],2) );
            
            
            % no laser =========================================================
            
            indPCD  = find( time >= 0 & time <= 0.20 );
            indBase = find( time >= -0.25 & time  <=0 );
            indDelay = find( time >=0.00 & time <=0.55 );
            
            PFC_muFR_PCD_C1{d}(i) = nanmean(spikeRateC1(indPCD));
            PFC_muFR_PCD_C2{d}(i) = nanmean(spikeRateC2(indPCD));
            
            PFC_muFR_CFD_C1{d}(i) = nanmean(spikeRateC1(indDelay));
            PFC_muFR_CFD_C2{d}(i) = nanmean(spikeRateC2(indDelay));
            
            PFC_muFR_BAS_C1{d}(i) = nanmean(spikeRateC1(indBase));
            PFC_muFR_BAS_C2{d}(i) = nanmean(spikeRateC2(indBase));
            
            PFC_CMI{d}(i)  =     CMI( PFC_muFR_PCD_C1{d}(i), PFC_muFR_PCD_C2{d}(i));
            PFC_CMI_BAS{d}(i)  = CMI( PFC_muFR_BAS_C1{d}(i), PFC_muFR_BAS_C2{d}(i));
            PFC_CMI_CFD{d}(i)  = CMI( PFC_muFR_CFD_C1{d}(i), PFC_muFR_CFD_C2{d}(i));
            
            PFC_CMI_C1_we{d}(i)  =  CMI( PFC_muFR_BAS_C1{d}(i), PFC_muFR_CFD_C1{d}(i));
            PFC_CMI_C2_we{d}(i)  =  CMI( PFC_muFR_BAS_C2{d}(i), PFC_muFR_CFD_C2{d}(i));
            
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            [spikeRateR1C2,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R1C2, bin, filtWidth, pre, post );
            
            PFC_muFR_PCD_R1C1{d}(i) = nanmean(spikeRateR1C1(indPCD));
            PFC_muFR_PCD_R1C2{d}(i) = nanmean(spikeRateR1C2(indPCD));
            
            PFC_muFR_CFD_R1C1{d}(i) = nanmean(spikeRateR1C1(indDelay));
            PFC_muFR_CFD_R1C2{d}(i) = nanmean(spikeRateR1C2(indDelay));
            
            PFC_muFR_BAS_R1C1{d}(i) = nanmean(spikeRateR1C1(indBase));
            PFC_muFR_BAS_R1C2{d}(i) = nanmean(spikeRateR1C2(indBase));
            
            
            [spikeRateR2C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            [spikeRateR2C2,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C2, bin, filtWidth, pre, post );
            
            PFC_muFR_PCD_R2C1{d}(i) = nanmean(spikeRateR2C1(indPCD));
            PFC_muFR_PCD_R2C2{d}(i) = nanmean(spikeRateR2C2(indPCD));
            
            PFC_muFR_CFD_R2C1{d}(i) = nanmean(spikeRateR2C1(indDelay));
            PFC_muFR_CFD_R2C2{d}(i) = nanmean(spikeRateR2C2(indDelay));
            
            PFC_muFR_BAS_R2C1{d}(i) = nanmean(spikeRateR2C1(indBase));
            PFC_muFR_BAS_R2C2{d}(i) = nanmean(spikeRateR2C2(indBase));
            
            
            PFC_Sel1{d}(i) =  PFC_muFR_CFD_R1C1{d}(i) -  PFC_muFR_CFD_R2C1{d}(i);
            PFC_Sel2{d}(i) =  PFC_muFR_CFD_R1C2{d}(i) -  PFC_muFR_CFD_R2C2{d}(i);
            
            RPFC{d}(i)   = sqrt( PFC_Sel1{d}(i).^2  +  PFC_Sel2{d}(i).^2 );
            ThetaPFC{d}(i) = atan( PFC_Sel1{d}(i)./ PFC_Sel2{d}(i) );
            
            
            % laser ==========================================================
            PFC_muFR_PCD_C1_las{d}(i) = nanmean(spikeRateC1Las(indPCD));
            PFC_muFR_PCD_C2_las{d}(i) = nanmean(spikeRateC2Las(indPCD));
            
            PFC_muFR_CFD_C1_las{d}(i) = nanmean(spikeRateC1Las(indDelay));
            PFC_muFR_CFD_C2_las{d}(i) = nanmean(spikeRateC2Las(indDelay));
            
            PFC_muFR_BAS_C1_las{d}(i) = nanmean(spikeRateC1Las(indBase));
            PFC_muFR_BAS_C2_las{d}(i) = nanmean(spikeRateC2Las(indBase));
            
            PFC_CMI_las{d}(i)  = CMI( PFC_muFR_PCD_C1_las{d}(i), PFC_muFR_PCD_C2_las{d}(i));
            PFC_CMI_BAS_las{d}(i)  = CMI( PFC_muFR_BAS_C1_las{d}(i), PFC_muFR_BAS_C2_las{d}(i));
            PFC_CMI_CFD_las{d}(i)  = CMI( PFC_muFR_CFD_C1_las{d}(i), PFC_muFR_CFD_C2_las{d}(i));
            
            PFC_CMI_C1_we_las{d}(i)  =  CMI( PFC_muFR_BAS_C1_las{d}(i), PFC_muFR_PCD_C1_las{d}(i));
            PFC_CMI_C2_we_las{d}(i)  =  CMI( PFC_muFR_BAS_C2_las{d}(i), PFC_muFR_PCD_C2_las{d}(i));
            
            
        [spikeRateR1C1Las,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1Las, bin, filtWidth, pre, post );
        [spikeRateR1C2Las,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R1C2Las, bin, filtWidth, pre, post );
        
        PFC_muFR_PCD_R1C1Las{d}(i) = nanmean(spikeRateR1C1Las(indPCD));
        PFC_muFR_PCD_R1C2Las{d}(i) = nanmean(spikeRateR1C2Las(indPCD));
        
        PFC_muFR_CFD_R1C1Las{d}(i) = nanmean(spikeRateR1C1Las(indDelay));
        PFC_muFR_CFD_R1C2Las{d}(i) = nanmean(spikeRateR1C2Las(indDelay));
        
        PFC_muFR_BAS_R1C1Las{d}(i) = nanmean(spikeRateR1C1Las(indBase));
        PFC_muFR_BAS_R1C2Las{d}(i) = nanmean(spikeRateR1C2Las(indBase));
        
        
        [spikeRateR2C1Las,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R2C1Las, bin, filtWidth, pre, post );
        [spikeRateR2C2Las,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C2Las, bin, filtWidth, pre, post );
        
        PFC_muFR_PCD_R2C1Las{d}(i) = nanmean(spikeRateR2C1Las(indPCD));
        PFC_muFR_PCD_R2C2Las{d}(i) = nanmean(spikeRateR2C2Las(indPCD));
        
        PFC_muFR_CFD_R2C1Las{d}(i) = nanmean(spikeRateR2C1Las(indDelay));
        PFC_muFR_CFD_R2C2Las{d}(i) = nanmean(spikeRateR2C2Las(indDelay));
        
        PFC_muFR_BAS_R2C1Las{d}(i) = nanmean(spikeRateR2C1Las(indBase));
        PFC_muFR_BAS_R2C2Las{d}(i) = nanmean(spikeRateR2C2Las(indBase));
        
        
        PFC_Sel1Las{d}(i) =  PFC_muFR_CFD_R1C1Las{d}(i) -  PFC_muFR_CFD_R2C1Las{d}(i);
        PFC_Sel2Las{d}(i) =  PFC_muFR_CFD_R1C2Las{d}(i) -  PFC_muFR_CFD_R2C2Las{d}(i);
        
        RPFCLas{d}(i)   = sqrt( PFC_Sel1Las{d}(i).^2  +  PFC_Sel2Las{d}(i).^2 );
        ThetaPFCLas{d}(i) = atan( PFC_Sel1Las{d}(i)./ PFC_Sel2Las{d}(i) );

        catch err
        end
    end;
end;

%%
for d = 1:length(ExptDates)
    
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} '/BPM33_catchTrials_Session.mat'] )
    
    %%
    
    [Spfc_mix, Smd_mix, ZVis, ZAud, first, ErrorFrac] = packageData_Laser(Z_C1, Smd, Spfc);
    FirstContext{d} = first;
    EF{d} = ErrorFrac;
    ct = 0;
    for i = 1:numel(Smd_mix)
        
        try
            
            clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateR1C1 spikeRateR2C1 spikeRateR1C2 spikeRateR2C2
            % correct nolaser =====================================================
            R1C1  = Smd_mix(i).SpikeTimes_R1_sound_NoLaser;
            R2C1  = Smd_mix(i).SpikeTimes_R2_sound_NoLaser;
            
            R1C2  = Smd_mix(i).SpikeTimes_R1_nosound_NoLaser;
            R2C2  = Smd_mix(i).SpikeTimes_R2_nosound_NoLaser;
            
            R1C1  =  R1C1( goodT(R1C1) );
            R1C2  =  R1C2( goodT(R1C2) );
            R2C1  =  R2C1( goodT(R2C1) );
            R2C2  =  R2C2( goodT(R2C2) );
            
            
            C1 = [R1C1, R2C1];
            C2 = [R1C2, R2C2];
            
            [spikeRateC1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
            [spikeRateC2,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
            
            ST = [C1, C2 ];
            [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(ST);
            base = mean_SC( find(bins_for_plotting < 0) );
            resp = mean_SC( find(bins_for_plotting > 0) );
            MDp{d}(i) = ranksum( base, resp );
            
            
            % correct laser =====================================================
            R1C1Las  = Smd_mix(i).SpikeTimes_R1_sound_Laser;
            R2C1Las  = Smd_mix(i).SpikeTimes_R2_sound_Laser;
            
            R1C2Las  = Smd_mix(i).SpikeTimes_R1_nosound_Laser;
            R2C2Las  = Smd_mix(i).SpikeTimes_R2_nosound_Laser;
            
            R1C1Las  =  R1C1Las;
            R1C2Las  =  R1C2Las;
            R2C1Las  =  R2C1Las;
            R2C2Las  =  R2C2Las;
            
            
            C1Las = [R1C1Las, R2C1Las];
            C2Las = [R1C2Las, R2C2Las];
            
            ST = [C1Las, C2Las ];
            [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(ST);
            base = mean_SC( find(bins_for_plotting < 0) );
            resp = mean_SC( find(bins_for_plotting > 0) );
            MDpLas{d}(i) = ranksum( base, resp );
            
            
            [spikeRateC1Las,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( C1Las, bin, filtWidth, pre, post );
            [spikeRateC2Las,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( C2Las, bin, filtWidth, pre, post );
            
            % no laser =========================================================
            
            indPCD  = find( time >= 0 & time <= 0.20 );
            indBase = find( time >= -0.25 & time  <=0 );
            indDelay = find( time >= 0.00 & time <=0.55 );
            
            MD_muFR_PCD_C1{d}(i) = nanmean(spikeRateC1(indPCD));
            MD_muFR_PCD_C2{d}(i) = nanmean(spikeRateC2(indPCD));
            
            MD_muFR_CFD_C1{d}(i) = nanmean(spikeRateC1(indDelay));
            MD_muFR_CFD_C2{d}(i) = nanmean(spikeRateC2(indDelay));
            
            MD_muFR_BAS_C1{d}(i) = nanmean(spikeRateC1(indBase));
            MD_muFR_BAS_C2{d}(i) = nanmean(spikeRateC2(indBase));
            
            MD_CMI{d}(i)  =     CMI( MD_muFR_PCD_C1{d}(i), MD_muFR_PCD_C2{d}(i));
            MD_CMI_BAS{d}(i)  = CMI( MD_muFR_BAS_C1{d}(i), MD_muFR_BAS_C2{d}(i));
            MD_CMI_CFD{d}(i)  = CMI( MD_muFR_CFD_C1{d}(i), MD_muFR_CFD_C2{d}(i));
            
            MD_CMI_C1_we{d}(i)  =  CMI( MD_muFR_BAS_C1{d}(i), MD_muFR_CFD_C1{d}(i));
            MD_CMI_C2_we{d}(i)  =  CMI( MD_muFR_BAS_C2{d}(i), MD_muFR_CFD_C2{d}(i));
            
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            [spikeRateR1C2,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R1C2, bin, filtWidth, pre, post );
            
            MD_muFR_PCD_R1C1{d}(i) = nanmean(spikeRateR1C1(indPCD));
            MD_muFR_PCD_R1C2{d}(i) = nanmean(spikeRateR1C2(indPCD));
            
            MD_muFR_CFD_R1C1{d}(i) = nanmean(spikeRateR1C1(indDelay));
            MD_muFR_CFD_R1C2{d}(i) = nanmean(spikeRateR1C2(indDelay));
            
            MD_muFR_BAS_R1C1{d}(i) = nanmean(spikeRateR1C1(indBase));
            MD_muFR_BAS_R1C2{d}(i) = nanmean(spikeRateR1C2(indBase));
            
            
            [spikeRateR2C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            [spikeRateR2C2,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C2, bin, filtWidth, pre, post );
            
            MD_muFR_PCD_R2C1{d}(i) = nanmean(spikeRateR2C1(indPCD));
            MD_muFR_PCD_R2C2{d}(i) = nanmean(spikeRateR2C2(indPCD));
            
            MD_muFR_CFD_R2C1{d}(i) = nanmean(spikeRateR2C1(indDelay));
            MD_muFR_CFD_R2C2{d}(i) = nanmean(spikeRateR2C2(indDelay));
            
            MD_muFR_BAS_R2C1{d}(i) = nanmean(spikeRateR2C1(indBase));
            MD_muFR_BAS_R2C2{d}(i) = nanmean(spikeRateR2C2(indBase));
            
            
            MD_Sel1{d}(i) =  MD_muFR_CFD_R1C1{d}(i) -  MD_muFR_CFD_R2C1{d}(i);
            MD_Sel2{d}(i) =  MD_muFR_CFD_R1C2{d}(i) -  MD_muFR_CFD_R2C2{d}(i);
            
            RMD{d}(i)   = sqrt( MD_Sel1{d}(i).^2  +  MD_Sel2{d}(i).^2 );
            ThetaMD{d}(i) = atan( MD_Sel1{d}(i)./ MD_Sel2{d}(i) );
            
            
            
            % laser ==========================================================
            MD_muFR_PCD_C1_las{d}(i) = nanmean(spikeRateC1Las(indPCD));
            MD_muFR_PCD_C2_las{d}(i) = nanmean(spikeRateC2Las(indPCD));
            
            MD_muFR_CFD_C1_las{d}(i) = nanmean(spikeRateC1Las(indDelay));
            MD_muFR_CFD_C2_las{d}(i) = nanmean(spikeRateC2Las(indDelay));
            
            MD_muFR_BAS_C1_las{d}(i) = nanmean(spikeRateC1Las(indBase));
            MD_muFR_BAS_C2_las{d}(i) = nanmean(spikeRateC2Las(indBase));
            
            MD_CMI_las{d}(i)  = CMI( MD_muFR_PCD_C1_las{d}(i), MD_muFR_PCD_C2_las{d}(i));
            MD_CMI_BAS_las{d}(i)  = CMI( MD_muFR_BAS_C1_las{d}(i), MD_muFR_BAS_C2_las{d}(i));
            MD_CMI_CFD_las{d}(i)  = CMI( MD_muFR_CFD_C1_las{d}(i), MD_muFR_CFD_C2_las{d}(i));
            
            MD_CMI_C1_we_las{d}(i)  =  CMI( MD_muFR_BAS_C1_las{d}(i), MD_muFR_PCD_C1_las{d}(i));
            MD_CMI_C2_we_las{d}(i)  =  CMI( MD_muFR_BAS_C2_las{d}(i), MD_muFR_PCD_C2_las{d}(i));
            
        [spikeRateR1C1Las,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1Las, bin, filtWidth, pre, post );
        [spikeRateR1C2Las,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R1C2Las, bin, filtWidth, pre, post );
        
        MD_muFR_PCD_R1C1Las{d}(i) = nanmean(spikeRateR1C1Las(indPCD));
        MD_muFR_PCD_R1C2Las{d}(i) = nanmean(spikeRateR1C2Las(indPCD));
        
        MD_muFR_CFD_R1C1Las{d}(i) = nanmean(spikeRateR1C1Las(indDelay));
        MD_muFR_CFD_R1C2Las{d}(i) = nanmean(spikeRateR1C2Las(indDelay));
        
        MD_muFR_BAS_R1C1Las{d}(i) = nanmean(spikeRateR1C1Las(indBase));
        MD_muFR_BAS_R1C2Las{d}(i) = nanmean(spikeRateR1C2Las(indBase));
        
        
        [spikeRateR2C1Las,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R2C1Las, bin, filtWidth, pre, post );
        [spikeRateR2C2Las,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C2Las, bin, filtWidth, pre, post );
        
        MD_muFR_PCD_R2C1Las{d}(i) = nanmean(spikeRateR2C1Las(indPCD));
        MD_muFR_PCD_R2C2Las{d}(i) = nanmean(spikeRateR2C2Las(indPCD));
        
        MD_muFR_CFD_R2C1Las{d}(i) = nanmean(spikeRateR2C1Las(indDelay));
        MD_muFR_CFD_R2C2Las{d}(i) = nanmean(spikeRateR2C2Las(indDelay));
        
        MD_muFR_BAS_R2C1Las{d}(i) = nanmean(spikeRateR2C1Las(indBase));
        MD_muFR_BAS_R2C2Las{d}(i) = nanmean(spikeRateR2C2Las(indBase));
        
        
        MD_Sel1Las{d}(i) =  MD_muFR_CFD_R1C1Las{d}(i) -  MD_muFR_CFD_R2C1Las{d}(i);
        MD_Sel2Las{d}(i) =  MD_muFR_CFD_R1C2Las{d}(i) -  MD_muFR_CFD_R2C2Las{d}(i);
        
        RMDLas{d}(i)   = sqrt( MD_Sel1Las{d}(i).^2  +  MD_Sel2Las{d}(i).^2 );
        ThetaMDLas{d}(i) = atan( MD_Sel1Las{d}(i)./ MD_Sel2Las{d}(i) );
            
        catch err
        end
    end;
    
    
    
end;

%%

colorMD = [0.49,0.18,0.56];
colorPFC = [0.5,0.5,0.5];
colorPFCLaser = [1, 0, 0];


PFCCMI_las = (cell2mat(PFC_CMI_CFD_las) );
PFCCMI     = (cell2mat( PFC_CMI_CFD)) ;

PFCCMI_las( PFCCMI_las == 1) = NaN;
PFCCMI_las( PFCCMI_las == -1) = NaN;

PFCCMI( PFCCMI == 1) = NaN;
PFCCMI( PFCCMI == -1) = NaN;


P = cell2mat( FracBlank );
PLas = cell2mat( FracBlankLaser );

gCells = find(P < 0.6);
gCellsLas  = find(PLas < 0.6);

figure(4); set(gcf,'color','w');

distributionPlot(PFCCMI(gCells)','globalNorm', 1, 'histOri','right','color','k','widthDiv',[2 2],'showMM',2)
distributionPlot(PFCCMI_las(gCellsLas)','globalNorm', 1,'histOri','left','color',colorPFC,'widthDiv',[2 1],'showMM',2)
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
ylim([-1,1]);



%%
MDCMI_las = (cell2mat(MD_CMI_CFD_las) );
MDCMI     = (cell2mat( MD_CMI_CFD));

MDCMI_las( MDCMI_las == 1 ) = NaN;
MDCMI_las( MDCMI_las == -1 ) = NaN;

MDCMI( MDCMI == 1) = NaN;
MDCMI( MDCMI == -1) = NaN;


P = cell2mat( MDp );
PLas = cell2mat( MDpLas );

gCells = find(P < 0.5);
gCellsLas  = find(PLas < 0.5);

figure(4); set(gcf,'color','w');

distributionPlot(MDCMI_las','globalNorm', 1, 'histOri','right','color',colorMD,'widthDiv',[2 2],'showMM',2)
distributionPlot(MDCMI','globalNorm', 1,'histOri','left','color',colorPFC,'widthDiv',[2 1],'showMM',2)
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
ylim([-1,1]);


%%


gCells = find(P < 0.05);


x = cell2mat( MD_Sel2 );
y = cell2mat( MD_Sel2Las );
ThetaC2 = atand( y./x );
RC2 = sqrt( x.^2 + y.^2 );


x = cell2mat( MD_Sel1 );
y = cell2mat( MD_Sel1Las );
ThetaC1 = atand( y./x );
RC1 = sqrt( x.^2 + y.^2 );

figure(100); set(gcf,'color','w');
 histogram( ThetaC2(RC2>1) , linspace(-90, 90, 20),'normalization','probability'); hold on;
 histogram( ThetaC1(RC1>1) , linspace(-90, 90, 20),'normalization','probability')
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')


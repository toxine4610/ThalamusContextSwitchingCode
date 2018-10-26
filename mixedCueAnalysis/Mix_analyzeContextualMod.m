
addpath ../Rajeev_Code
% Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/BPM1_3/ContextSwitchHalo/UniLatMD-Low/';

%%
CMI   = @(X,Y) (X-Y)./(X+Y);
goodT =  @(X) find(~cellfun(@isempty,X));
bin = 0.01;
filtWidth = 0.05;
pre =  0.2;
post = 1.8;

%%
%  ExptDates = {'2017-12-15','2017-12-16','2017-12-17','2017-12-18','2017-12-20','2017-12-22','2017-12-23','2017-12-24','2017-12-26'};
%  ExptDates = {'2018-03-13','2018-03-14','2018-03-15','2018-03-16','2018-03-17', '2018-03-18','2018-03-20','2018-03-22'};

% ExptDates = {'2017-12-15'} %,'2017-12-16','2017-12-17','2017-12-18','2017-12-20','2017-12-22','2017-12-23','2017-12-24','2017-12-26'};

Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/BPM3_3/CatchTrials/';
ExptDates = {'2018-04-08', '2018-04-09', '2018-04-10'};


%%
for d = 1:length(ExptDates)
     
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} filesep 'BPM33_catchTrials_Session.mat'] )
    
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
        try
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
        
        [spikeRateC1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
        [spikeRateC2,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
        
        PFCFRTS_C1{d}(i,:) = spikeRateC1;
        PFCFRTS_C2{d}(i,:) = spikeRateC2;
        
%         [C1_weightedCC, C1_weightedCC_scaled{d}(i,:)] = computeRawReliability(C1);
%         [C2_weightedCC, C2_weightedCC_scaled{d}(i,:)] = computeRawReliability(C2);
%         
        indPCD  = find( time >= 0 & time <= 0.25 );
        indBase = find( time >= -0.25 & time  <=0 );
        indDelay = find( time >=0.35 & time <=0.85 );
        
        PFC_muFR_PCD_C1{d}(i) = nanmean(spikeRateC1(indPCD));
        PFC_muFR_PCD_C2{d}(i) = nanmean(spikeRateC2(indPCD));
        
        PFC_muFR_CFD_C1{d}(i) = nanmean(spikeRateC1(indDelay));
        PFC_muFR_CFD_C2{d}(i) = nanmean(spikeRateC2(indDelay));
        
        PFC_muFR_BAS_C1{d}(i) = nanmean(spikeRateC1(indBase));
        PFC_muFR_BAS_C2{d}(i) = nanmean(spikeRateC2(indBase));
        
        PFC_CMI{d}(i)  =     CMI( PFC_muFR_PCD_C1{d}(i), PFC_muFR_PCD_C2{d}(i));
        PFC_CMI_BAS{d}(i)  = CMI( PFC_muFR_BAS_C1{d}(i), PFC_muFR_BAS_C2{d}(i));
        PFC_CMI_CFD{d}(i)  = CMI( PFC_muFR_CFD_C1{d}(i), PFC_muFR_CFD_C2{d}(i));
        
        PFC_CMI_C1_we{d}(i)  =  CMI( PFC_muFR_PCD_C1{d}(i), PFC_muFR_CFD_C1{d}(i));
        PFC_CMI_C2_we{d}(i)  =  CMI( PFC_muFR_PCD_C2{d}(i), PFC_muFR_CFD_C2{d}(i));
        %
        if size(R1C1,2) ~= 0
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            SM11 =  spikeMatR1C1';
          
        else
            spikeRateR1C1 = NaN(1,201);
            SM11 = NaN(22, 201) ;
             
        end;
        
        if size(R1C2,2) ~= 0
          [spikeRateR1C2,errorBounds,spikeMatR1C2,time] = estimateSpikeRates( R1C2, bin, filtWidth, pre, post );
          SM12 =  spikeMatR1C2';
        else
          spikeRateR1C2 = NaN(1,201);
          SM12 = NaN(22,201);
        end;
            
        PFC_muFR_PCD_R1C1{d}(i) = nanmean(spikeRateR1C1(indPCD));
        PFC_muFR_PCD_R1C2{d}(i) = nanmean(spikeRateR1C2(indPCD));
        
        PFC_muFR_CFD_R1C1{d}(i) = nanmean(spikeRateR1C1(indDelay));
        PFC_muFR_CFD_R1C2{d}(i) = nanmean(spikeRateR1C2(indDelay));
        
        PFC_muFR_BAS_R1C1{d}(i) = nanmean(spikeRateR1C1(indBase));
        PFC_muFR_BAS_R1C2{d}(i) = nanmean(spikeRateR1C2(indBase));
        
        
        if size(R2C1,2) ~= 0
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            SM21 =  spikeMatR2C1';
        else
            spikeRateR2C1 = NaN(1,201);
            SM21 = NaN(22, 201);
        end;
        
        if size(R2C2,2) ~= 0
          [spikeRateR2C2,errorBounds,spikeMatR2C2,time] = estimateSpikeRates( R2C2, bin, filtWidth, pre, post );
          SM22 =  spikeMatR2C2';
        else
            spikeRateR2C2 = NaN(1,201);
            SM22 = NaN(22,201);
        end;
        if length(spikeRateR2C1) == 1
            spikeRateR2C1 = NaN(1,length(time));
        end;
        if length(spikeRateR2C2) == 1
            spikeRateR2C2 = NaN(1,length(time));
        end;
        
        PFC_muFR_PCD_R2C1{d}(i) = nanmean(spikeRateR2C1(indPCD));
        PFC_muFR_PCD_R2C2{d}(i) = nanmean(spikeRateR2C2(indPCD));
        
        PFC_muFR_CFD_R2C1{d}(i) = nanmean(spikeRateR2C1(indDelay));
        PFC_muFR_CFD_R2C2{d}(i) = nanmean(spikeRateR2C2(indDelay));
        
        PFC_muFR_BAS_R2C1{d}(i) = nanmean(spikeRateR2C1(indBase));
        PFC_muFR_BAS_R2C2{d}(i) = nanmean(spikeRateR2C2(indBase));
        
        PFC_Sel1{d}(i) =  ( PFC_muFR_CFD_R1C1{d}(i) -  PFC_muFR_CFD_R2C1{d}(i));
        PFC_Sel2{d}(i) =  (PFC_muFR_CFD_R1C2{d}(i) -  PFC_muFR_PCD_R2C2{d}(i));
        
        RPFC{d}(i)   = sqrt( PFC_Sel1{d}(i).^2  +  PFC_Sel2{d}(i).^2 );
        ThetaPFC{d}(i) = atan( PFC_Sel1{d}(i)./ PFC_Sel2{d}(i) );
        
        
        
        PFC_dp_C1{d}(i,:) = (nanmean( SM11 ) - nanmean( SM21 ))./ sqrt( 0.5.*( var( SM21,[],1) + var( SM11,[], 1) ) );
        PFC_dp_C2{d}(i,:) = (nanmean( SM12 ) - nanmean( SM22 ))./ sqrt( 0.5.*( var( SM12,[],1) + var( SM22,[], 1) ) );

        
        
        
        %
        %
        %         PFC_CMI_R2{d}(i)  = CMI( PFC_muFR_PCD_R2C1{d}(i), PFC_muFR_PCD_R2C2{d}(i));
        %
        %         PFC_CMI_BAS_R1{d}(i)  = CMI( PFC_muFR_BAS_R1C1{d}(i), PFC_muFR_BAS_R1C2{d}(i));
        %         PFC_CMI_BAS_R2{d}(i)  = CMI( PFC_muFR_BAS_R2C1{d}(i), PFC_muFR_BAS_R2C2{d}(i));
        %
        %         PFC_CMI_CFD_R1{d}(i)  = CMI( PFC_muFR_CFD_R1C1{d}(i), PFC_muFR_CFD_R1C2{d}(i));
        %         PFC_CMI_CFD_R2{d}(i)  = CMI( PFC_muFR_CFD_R2C1{d}(i), PFC_muFR_CFD_R2C2{d}(i));
        %
        catch err
        end
    end;
    
    %%
    for i = setdiff( 1:numel(Smd_mix), [] )
        
        clear C1 C2 R1C1 R2C1 spikeRateR1C1 spikeRateR2C1 RC2 R2C2 spikeRateR1C2 spikeRateR2C2
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
        
        
        [spikeRateC1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
        [spikeRateC2,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
        
        
        MDFRTS_C1{d}(i,:) = spikeRateC1;
        MDFRTS_C2{d}(i,:) = spikeRateC2;
        
        
        indPCD  = find( time >= 0 & time <= 0.25 );
        indBase = find( time >= -0.25 & time  <=0 );
        indDelay = find( time >=0.35 & time <=0.85 );
        
        MD_muFR_PCD_C1{d}(i) = nanmean(spikeRateC1(indPCD));
        MD_muFR_PCD_C2{d}(i) = nanmean(spikeRateC2(indPCD));
        
        MD_muFR_CFD_C1{d}(i) = nanmean(spikeRateC1(indDelay));
        MD_muFR_CFD_C2{d}(i) = nanmean(spikeRateC2(indDelay));
        
        MD_muFR_BAS_C1{d}(i) = nanmean(spikeRateC1(indBase));
        MD_muFR_BAS_C2{d}(i) = nanmean(spikeRateC2(indBase));
        
        MD_CMI{d}(i)      = CMI( MD_muFR_PCD_C1{d}(i), MD_muFR_PCD_C2{d}(i));
        MD_CMI_BAS{d}(i)  = CMI( MD_muFR_BAS_C1{d}(i), MD_muFR_BAS_C2{d}(i));
        MD_CMI_CFD{d}(i)  = CMI( MD_muFR_CFD_C1{d}(i), MD_muFR_CFD_C2{d}(i));
        
        MD_CMI_C1_we{d}(i)  = CMI( MD_muFR_PCD_C1{d}(i), MD_muFR_CFD_C1{d}(i));
        MD_CMI_C2_we{d}(i)  = CMI( MD_muFR_PCD_C2{d}(i), MD_muFR_CFD_C2{d}(i));
        
        %
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
        %
        %
        %         MD_CMI_R1{d}(i)  = CMI( MD_muFR_PCD_R1C1{d}(i), MD_muFR_PCD_R1C2{d}(i));
        %         MD_CMI_R2{d}(i)  = CMI( MD_muFR_PCD_R2C1{d}(i), MD_muFR_PCD_R2C2{d}(i));
        %
        %         MD_CMI_BAS_R1{d}(i)  = CMI( MD_muFR_BAS_R1C1{d}(i), MD_muFR_BAS_R1C2{d}(i));
        %         MD_CMI_BAS_R2{d}(i)  = CMI( MD_muFR_BAS_R2C1{d}(i), MD_muFR_BAS_R2C2{d}(i));
        %
        %         MD_CMI_CFD_R1{d}(i)  = CMI( MD_muFR_CFD_R1C1{d}(i), MD_muFR_CFD_R1C2{d}(i));
        %         MD_CMI_CFD_R2{d}(i)  = CMI( MD_muFR_CFD_R2C1{d}(i), MD_muFR_CFD_R2C2{d}(i));
        
    end;
end;


%%

colorMD = [0.49,0.18,0.56];
colorPFC = [0.5,0.5,0.5];

MDCMI  = (cell2mat(MD_CMI) );
PFCCMI = (cell2mat( PFC_CMI)) ;

% MD_BAS_CMI  = nanmean( [MD_CMI_BAS_R2{1}; MD_CMI_BAS_R1{1}] );
% PFC_BAS_CMI = nanmean( [PFC_CMI_BAS_R2{1}; PFC_CMI_BAS_R1{1}] );
%
% MD_CFD_CMI  = nanmean( [MD_CMI_CFD_R2{1}; MD_CMI_CFD_R1{1}] );
% PFC_CFD_CMI = nanmean( [PFC_CMI_CFD_R2{1}; PFC_CMI_CFD_R1{1}] );

figure(4); set(gcf,'color','w');
% subplot(1,3,1)
% distributionPlot(PFC_BAS_CMI','globalNorm', 1, 'histOri','right','color',colorPFC,'widthDiv',[2 2],'showMM',2)
% distributionPlot(MD_BAS_CMI','globalNorm', 1,'histOri','left','color',colorMD,'widthDiv',[2 1],'showMM',2)
% axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
% ylim([-1,1]);

% subplot(1,3,2)
distributionPlot(PFCCMI','globalNorm', 1, 'histOri','right','color',colorPFC,'widthDiv',[2 2],'showMM',2)
distributionPlot(MDCMI','globalNorm', 1,'histOri','left','color',colorMD,'widthDiv',[2 1],'showMM',2)
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
ylim([-1,1]);


% figure(100);
% datacell = {MDCMI , PFCCMI };
% plotSpread(datacell, 'categoryMarkers',{'o','o'},'categoryColors',{colorMD,colorPFC});

% subplot(1,3,3)
% distributionPlot(PFC_CFD_CMI','globalNorm', 1, 'histOri','right','color',colorPFC,'widthDiv',[2 2],'showMM',2)
% distributionPlot(MD_CFD_CMI','globalNorm', 1,'histOri','left','color',colorMD,'widthDiv',[2 1],'showMM',2)
% axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
% ylim([-1,1]);



%%
figure(4); set(gcf,'color','w');
plot( MD_muFR_CFD_R1C1{1}, MD_muFR_CFD_R1C2{1},'o','markersize',12,'markerfacecolor',colorMD,'markeredgecolor','w');
hold on;
plot( PFC_muFR_CFD_R2C1{1}, PFC_muFR_CFD_R2C2{1},'o','markersize',12,'markerfacecolor',colorPFC,'markeredgecolor','w');
hold on;

plot( MD_muFR_PCD_R1C1{1}, MD_muFR_PCD_R1C2{1},'o','markersize',12,'markerfacecolor',colorMD,'markeredgecolor','w');
hold on;
plot( PFC_muFR_PCD_R2C1{1}, PFC_muFR_PCD_R2C2{1},'o','markersize',12,'markerfacecolor',colorPFC,'markeredgecolor','w');
hold on;

plot( linspace(0,12,2000), linspace(0,12,2000),'color','k');
axis square; box off; set(gca,'tickdir','out','fontsize',16);
set(gca,'xscale','log','yscale','log')


%%

figure(5); set(gcf,'color','w');

MDS1  = abs(cell2mat(PFC_Sel1));
MDS2  = abs(cell2mat(PFC_Sel2));
ind1  = find( MDS1 >= 0.5 | MDS1 <= -0.5);
ind2  = find( MDS2 >= 0.5 | MDS2 <= -0.5 );
ind   = intersect(ind1, ind2);

plot( MDS1(ind), MDS2(ind),'o','markersize',12,'markerfacecolor',colorMD,'markeredgecolor','w');
hold on;
axis square; box off; set(gca,'tickdir','out','fontsize',16)

SelAngle = atand( MDS2(ind)./ MDS1(ind) );
SelStrength = sqrt( MDS1(ind).^2  +  MDS2(ind).^2 );

figure(7); set(gcf,'color','w');
histogram(SelAngle, linspace(-90,90,20),'normalization','probability'); hold on;
axis normal; box off; set(gca,'tickdir','out','fontsize',16);



%%

figure(6);
for i = 1:d
    plot( PerfDiff{i}, nanmedian(PFC_CMI{i}), 'o','markersize',12,'markerfacecolor',colorPFC,'markeredgecolor','w');
    hold on;
    plot( PerfDiff{i}, nanmedian(MD_CMI{i}), 'o','markersize',12,'markerfacecolor',colorMD,'markeredgecolor','w');
end;
    


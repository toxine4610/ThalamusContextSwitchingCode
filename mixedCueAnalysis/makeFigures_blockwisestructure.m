FigSaveDir = [D '/'];
mkdir(FigSaveDir);
ANA = 0;

laserBlack = [0.5, 0.5, 0.5];
laserRed   = [219./255, 112./255, 147./255];

%%

indCorrVis = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 1);
indCorrAud = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 2);

ZVis = Z_C1(indCorrVis, :);
ZAud = Z_C1(indCorrAud, :);

indVisNoSound = find( ZVis(:,9) == 3 & ZVis(:,10) == 0);
indAudNoSound = find( ZAud(:, 9) == 3 & ZAud(:,10) == 0);

indVisSound = find( ZVis(:,9) == 0 & ZVis(:,10) == 0);
indAudSound = find( ZAud(:, 9) == 0 & ZAud(:,10) == 0);

indVisNoSound_Laser = find( ZVis(:,9) == 3 & ZVis(:,10) == 1);
indAudNoSound_Laser = find( ZAud(:, 9) == 3 & ZAud(:,10) == 1);

indVisSound_Laser = find( ZVis(:,9) == 0 & ZVis(:,10) == 1);
indAudSound_Laser = find( ZAud(:, 9) == 0 & ZAud(:,10) == 1);


indBadVis = find(Z_C1(:,1) == 0 & Z_C1(:,2) == 1);
indBadAud = find(Z_C1(:,1) == 0 & Z_C1(:,2) == 2);

ZVisBad = Z_C1(indBadVis, :);
ZAudBad = Z_C1(indBadAud, :);

indVisNoSoundBad = find( ZVisBad(:,9) == 3 & ZVisBad(:,10)==0);
indAudNoSoundBad = find( ZAudBad(:, 9) == 3 & ZAudBad(:,10)==0);

indVisSoundBad = find( ZVisBad(:,9) == 0 & ZVisBad(:,10)==0);
indAudSoundBad = find( ZAudBad(:, 9) == 0 & ZAudBad(:,10)==0);

indVisNoSoundBad_Laser = find( ZVisBad(:,9) == 3 & ZVisBad(:,10)==1);
indAudNoSoundBad_Laser = find( ZAudBad(:, 9) == 3 & ZAudBad(:,10)==1);

indVisSoundBad_Laser = find( ZVisBad(:,9) == 0 & ZVisBad(:,10)==1);
indAudSoundBad_Laser = find( ZAudBad(:, 9) == 0 & ZAudBad(:,10)==1);

%%
indAudBlockLaser = find(Z_C1(:,10)==1 & Z_C1(:,9)==0);
indVisBlockLaser = find(Z_C1(:,10)==1 & Z_C1(:,9)==3);

indAudBlockNoLaser = find(Z_C1(:,10)==0 & Z_C1(:,9)==0);
indVisBlockNoLaser = find(Z_C1(:,10)==0 & Z_C1(:,9)==3);

foo  = Z_C1(indAudBlockLaser,:);
foo2 = Z_C1(indAudBlockNoLaser,:);
foo3 = Z_C1(indVisBlockNoLaser,:);
foo4  = Z_C1(indVisBlockLaser,:);

ErrorFrac_Aud_Block_Laser = length( find(foo(:,1)==0) )./size(foo,1);
ErrorFrac_Aud_Block_NoLaser = length( find(foo2(:,1)==0) )./size(foo2,1);
ErrorFrac_Vis_Block_NoLaser = length( find(foo3(:,1)==0) )./size(foo3,1);
ErrorFrac_Vis_Block_Laser = length( find(foo4(:,1)==0) )./size(foo4,1);

figure(200); set(gcf,'color','w')
b=bar(1, ErrorFrac_Aud_Block_NoLaser); hold on;
set(b(1),'facecolor','k');

b=bar(2, ErrorFrac_Aud_Block_Laser); hold on;
set(b(1),'facecolor',laserBlack);

b=bar(3, ErrorFrac_Vis_Block_NoLaser); hold on;
set(b(1),'facecolor','r');

b=bar(4, ErrorFrac_Vis_Block_Laser); hold on;
set(b(1),'facecolor',laserRed);

ylabel('Overall Error fraction')
axis square; box off; set(gca,'tickdir','out','fontsize',16);

ylim([0,0.6]);
line([0, 4],[0.5,0.5]);


%% Analyze performance
if ANA == 1
    Performance = Z_C1(:,8);
    trialType = Z_C1(:,9);
    d = diff( trialType );
    Start =  find(d==3)+1; End = find( d==-3 );
    V = [Start, End];
    
    for i = 1:size(V,1)
        s = V(i,1); e = V(i,2);
        nC(i) = length( find( Z_C1( s:e, 1) == 1));
        nW(i) = length( find( Z_C1( s:e, 1) == 0));
        F(i)     = nW(i)./(nC(i)+nW(i));
        nT(i) = (nC(i)+nW(i));
    end;
    
    for i = 1:size(V,1)
        if i == 1
            s = 1;  e = V(i,1)-1;
            nC(i) = length( find( Z_C1( s:e, 1) == 1));
            nW(i) = length( find( Z_C1( s:e, 1) == 0));
            Fs(i)     = nW(i)./(nC(i)+nW(i));
            nTs(i) = nC(i)+nW(i);
        elseif i >= 1
            s = V(i-1,2)+1; e = V(i,1)-1;
            nC(i) = length( find( Z_C1( s:e, 1) == 1));
            nW(i) = length( find( Z_C1( s:e, 1) == 0));
            Fs(i)     = nW(i)./(nC(i)+nW(i));
            nTs(i) = (nC(i)+nW(i));
        end;
    end;
    
    figure(400); set(gcf,'color','w');
    plot( [1:size(V,1)], F, '-o', 'markersize',12,'markerfacecolor','r','markeredgecolor','r','color','r'); hold on;
    plot( 1:size(V,1), Fs, '-o', 'markersize',12,'markerfacecolor','k','markeredgecolor','k','color','k');
    xlim([0, size(V,1)+1]);
    line( [0, size(V,1)+1], [0.5, 0.5],'color','k' );
    ylabel('Block-wise error fraction'); xlabel('Block Number');
    axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
    
    figure(500); set(gcf,'color','w');
    plot( 1:size(V,1), nT(1:end), '-o', 'markersize',12,'markerfacecolor','r','markeredgecolor','r','color','r'); hold on;
    plot( 1:size(V,1), nTs(1:end), '-o', 'markersize',12,'markerfacecolor','k','markeredgecolor','k','color','k');
    xlim([0, size(V,1)+1]);
    ylabel('Number of Trials to Complete Block'); xlabel('Block Number');
    axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
    
end;
%%

for i = 1:numel(Smd)
    
    Smd_mix(i).SpikeTimes_R1_sound = Smd(i).SpikeTimes_R1C1(indVisSound);
    Smd_mix(i).SpikeTimes_R2_sound = Smd(i).SpikeTimes_R2C1(indAudSound);
    
    Smd_mix(i).SpikeTimes_R1_nosound = Smd(i).SpikeTimes_R1C1(indVisNoSound);
    Smd_mix(i).SpikeTimes_R2_nosound = Smd(i).SpikeTimes_R2C1(indAudNoSound);
    
    
    Smd_mix(i).SpikeTimes_R1_sound_laser = Smd(i).SpikeTimes_R1C1(indVisSound_Laser);
    Smd_mix(i).SpikeTimes_R2_sound_laser = Smd(i).SpikeTimes_R2C1(indAudSound_Laser);
    
    Smd_mix(i).SpikeTimes_R1_nosound_laser = Smd(i).SpikeTimes_R1C1(indVisNoSound_Laser);
    Smd_mix(i).SpikeTimes_R2_nosound_laser = Smd(i).SpikeTimes_R2C1(indAudNoSound_Laser);
    
    
    Smd_mix(i).SpikeTimes_R1_sound_inc = Smd(i).SpikeTimes_R1C1_inc(indVisSoundBad);
    Smd_mix(i).SpikeTimes_R2_sound_inc = Smd(i).SpikeTimes_R2C1_inc(indAudSoundBad);
    
    Smd_mix(i).SpikeTimes_R1_nosound_inc = Smd(i).SpikeTimes_R1C1_inc(indVisNoSoundBad);
    Smd_mix(i).SpikeTimes_R2_nosound_inc = Smd(i).SpikeTimes_R2C1_inc(indAudNoSoundBad);
    
    Smd_mix(i).SpikeTimes_R1_sound_inc_laser = Smd(i).SpikeTimes_R1C1_inc(indVisSoundBad_Laser);
    Smd_mix(i).SpikeTimes_R2_sound_inc_laser = Smd(i).SpikeTimes_R2C1_inc(indAudSoundBad_Laser);
    
    Smd_mix(i).SpikeTimes_R1_nosound_inc_laser = Smd(i).SpikeTimes_R1C1_inc(indVisNoSoundBad_Laser);
    Smd_mix(i).SpikeTimes_R2_nosound_inc_laser = Smd(i).SpikeTimes_R2C1_inc(indAudNoSoundBad_Laser);
    
end;

for i = 1:numel(Spfc)
    
    Spfc_mix(i).SpikeTimes_R1_sound = Spfc(i).SpikeTimes_R1C1(indVisSound);
    Spfc_mix(i).SpikeTimes_R2_sound = Spfc(i).SpikeTimes_R2C1(indAudSound);
    
    Spfc_mix(i).SpikeTimes_R1_nosound = Spfc(i).SpikeTimes_R1C1(indVisNoSound);
    Spfc_mix(i).SpikeTimes_R2_nosound = Spfc(i).SpikeTimes_R2C1(indAudNoSound);
    
    Spfc_mix(i).SpikeTimes_R1_sound_laser = Spfc(i).SpikeTimes_R1C1(indVisSound_Laser);
    Spfc_mix(i).SpikeTimes_R2_sound_laser = Spfc(i).SpikeTimes_R2C1(indAudSound_Laser);
    
    Spfc_mix(i).SpikeTimes_R1_nosound_laser = Spfc(i).SpikeTimes_R1C1(indVisNoSound_Laser);
    Spfc_mix(i).SpikeTimes_R2_nosound_laser = Spfc(i).SpikeTimes_R2C1(indAudNoSound_Laser);
    
    
    Spfc_mix(i).SpikeTimes_R1_sound_inc = Spfc(i).SpikeTimes_R1C1_inc(indVisSoundBad);
    Spfc_mix(i).SpikeTimes_R2_sound_inc = Spfc(i).SpikeTimes_R2C1_inc(indAudSoundBad);
    
    Spfc_mix(i).SpikeTimes_R1_nosound_inc = Spfc(i).SpikeTimes_R1C1_inc(indVisNoSoundBad);
    Spfc_mix(i).SpikeTimes_R2_nosound_inc = Spfc(i).SpikeTimes_R2C1_inc(indAudNoSoundBad);
    
    
    Spfc_mix(i).SpikeTimes_R1_sound_inc_laser = Spfc(i).SpikeTimes_R1C1_inc(indVisSoundBad_Laser);
    Spfc_mix(i).SpikeTimes_R2_sound_inc_laser = Spfc(i).SpikeTimes_R2C1_inc(indAudSoundBad_Laser);
    
    Spfc_mix(i).SpikeTimes_R1_nosound_inc_laser = Spfc(i).SpikeTimes_R1C1_inc(indVisNoSoundBad_Laser);
    Spfc_mix(i).SpikeTimes_R2_nosound_inc_laser = Spfc(i).SpikeTimes_R2C1_inc(indAudNoSoundBad_Laser);
end;

goodT =  @(X) find(~cellfun(@isempty,X));

delayEnd = 1.2 ; preCueEnd = 0.2;


range = [ -0.3, 0.7 ];
bin = 0.0010;
filtWidth = 0.08;

xplotLim = [-0.25, 1.25];
yplotLim = [0,20];

%%
for i = 1:numel(Smd_mix)
    
    try
        
        PSTH_Raster_plot =  figure(i);
        PSTH_Raster_plot.Renderer='Painters';
        set(PSTH_Raster_plot, 'position', [63  39  1871  1318])
        set(PSTH_Raster_plot,'color','w');
        set(0,'DefaultAxesFontSize',16)
        set(0,'DefaultLineLineWidth',1)
        hold on;
        
        R1Sound = Smd_mix(i).SpikeTimes_R1_sound( goodT(Smd_mix(i).SpikeTimes_R1_sound));
        R2Sound = Smd_mix(i).SpikeTimes_R2_sound( goodT(Smd_mix(i).SpikeTimes_R2_sound));
        R1NoSound = Smd_mix(i).SpikeTimes_R1_nosound( goodT(Smd_mix(i).SpikeTimes_R1_nosound));
        R2NoSound = Smd_mix(i).SpikeTimes_R2_nosound( goodT(Smd_mix(i).SpikeTimes_R2_nosound));
        
        R1Sound_Laser = [ Smd_mix(i).SpikeTimes_R1_sound_laser, Smd_mix(i).SpikeTimes_R1_sound_inc_laser];
        R2Sound_Laser = [ Smd_mix(i).SpikeTimes_R2_sound_laser, Smd_mix(i).SpikeTimes_R2_sound_inc_laser];
        R1NoSound_Laser = [ Smd_mix(i).SpikeTimes_R1_nosound_laser, Smd_mix(i).SpikeTimes_R1_nosound_inc_laser];
        R2NoSound_Laser = [ Smd_mix(i).SpikeTimes_R2_nosound_laser, Smd_mix(i).SpikeTimes_R2_nosound_inc_laser];
        
        subplot(3,4,1);
        plotRaster( R1Sound );
        xlim(xplotLim);
        title('R1 - Aud','color','k');
        ylim( yplotLim );
        
        subplot(3,4,2);
        plotRaster( R2Sound );
        xlim(xplotLim);
        ylim( yplotLim );
        title('R2 - Aud','color','k');
        
        subplot(3,4,3);
        plotRaster( R1NoSound );
        xlim(xplotLim);
        ylim( yplotLim );
        title('R1 - Vis','color','r');
        
        subplot(3,4,4);
        plotRaster( R2NoSound );
        xlim(xplotLim);
        ylim( yplotLim );
        title('R2 - Vis','color','r');
        
        if size(R1Sound_Laser,2)~= 0
            subplot(3,4,5);
            plotRaster( R1Sound_Laser );
            xlim(xplotLim);
            title('R1 - Aud (Laser)','color',laserBlack);
            ylim( yplotLim );
        end;
        
        if size(R2Sound_Laser,2)~= 0
            subplot(3,4,6);
            plotRaster( R2Sound_Laser );
            xlim(xplotLim);
            ylim( yplotLim );
            title('R2 - Aud (Laser)','color',laserBlack);
        end;
        
        if size(R1NoSound_Laser,2)~= 0
            subplot(3,4,7);
            plotRaster( R1NoSound_Laser );
            xlim(xplotLim);
            ylim( yplotLim );
            title('R1 - Vis (Laser)','color',laserRed);
        end
        
        if size(R2NoSound_Laser)~=0
            subplot(3,4,8);
            plotRaster( R2NoSound_Laser );
            xlim(xplotLim);
            ylim( yplotLim );
            title('R2 - Vis (Laser)','color',laserRed);
        end
        
        % R1 - Aud ============================================================
        [R1SoundZ,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  R1Sound, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
        
        %         [R1SoundZ, spikeMat, time1, indBase] = makeSpikeRates(R1Sound, range , bin, filtWidth);
        
        subplot(3,4,[9,10]);
        plot( time1, R1SoundZ,'color','k','linewidth',3); hold on;
        box off; set(gca,'tickdir','out','fontsize',16);
        xlim(xplotLim);
        yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        title('Rule 1');
        
        % R1 - vis ============================================================
        [R1NoSoundZ,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  R1NoSound, 0.02, [-0.2 0.0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
        %         [R1NoSoundZ, spikeMat, time1, indBase] = makeSpikeRates(R1NoSound, range , bin, filtWidth);
        
        
        subplot(3,4,[9,10]);
        plot( time1, R1NoSoundZ,'color','r','linewidth',3);
        box off; set(gca,'tickdir','out','fontsize',16);
        xlim(xplotLim);     yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        
        if size(R1Sound_Laser,2)~= 0
            % R1 - Aud laser ============================================================
            %     [R1Sound_laser_Z,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  R1Sound_Laser, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
            [R1Sound_laser_Z, spikeMat, time1, indBase] = makeSpikeRates(R1Sound_Laser, range , bin, filtWidth);
            
            
            subplot(3,4,[9,10]);
            plot( time1, R1Sound_laser_Z,'color',laserBlack,'linewidth',3); hold on;
            box off; set(gca,'tickdir','out','fontsize',16);
            xlim(xplotLim);
            yLims = get(gca, 'YLim');
            line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        end
        
        if size(R1NoSound_Laser,2)~= 0
            % R1 - vis laser ============================================================
            %     [R1NoSound_laser_Z,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  R1NoSound_Laser, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
            [R1NoSound_laser_Z, spikeMat, time1, indBase] = makeSpikeRates(R1NoSound_Laser, range , bin, filtWidth);
            
            subplot(3,4,[9,10]);
            plot( time1, R1NoSound_laser_Z,'color',laserRed,'linewidth',3);
            box off; set(gca,'tickdir','out','fontsize',16);
            xlim(xplotLim);
            yLims = get(gca, 'YLim');
            line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        end;
        
        
        
        % R2 - Aud ============================================================
        %             [R2SoundZ,spikeMatFilt2,baseStd,baseMean,time2] = spikeZscoreSigTest( R2Sound, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
        [R2SoundZ, spikeMat, time1, indBase] = makeSpikeRates(R2Sound, range , bin, filtWidth);
        
        
        subplot(3,4,[11,12]);
        plot( time1, R2SoundZ,'color','k','linewidth',3);hold on;
        box off; set(gca,'tickdir','out','fontsize',16);
        xlim(xplotLim);
        yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        title('Rule 2');
        
        % R2 - Vis ============================================================
        %              [R2NoSoundZ,spikeMatFilt2,baseStd,baseMean,time2] = spikeZscoreSigTest( R2NoSound, 0.02,[-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
        [R2NoSoundZ, spikeMat, time1, indBase] = makeSpikeRates(R2NoSound, range , bin, filtWidth);
        
        
        subplot(3,4,[11,12]);
        plot( time1, R2NoSoundZ,'color','r','linewidth',3);
        box off; set(gca,'tickdir','out','fontsize',16);
        xlim(xplotLim);
        yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        
        
        if size(R2Sound_Laser,2)~= 0
            % R2 - Aud - laser ============================================================
            %     [R2Sound_laser_Z,spikeMatFilt2,baseStd,baseMean,time2] = spikeZscoreSigTest( R2Sound_Laser, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
            [R2Sound_laser_Z, spikeMat, time1, indBase] = makeSpikeRates(R2Sound_Laser, range , bin, filtWidth);
            
            subplot(3,4,[11,12]);
            plot( time1, R2Sound_laser_Z,'color',laserBlack,'linewidth',3);hold on;
            box off; set(gca,'tickdir','out','fontsize',16);
            xlim(xplotLim);
            yLims = get(gca, 'YLim');
            line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
            
        end;
        
        if size(R2NoSound_Laser,2)~= 0
            % R2 - Vis - laser ============================================================
            %     [R2NoSound_laser_Z,spikeMatFilt2,baseStd,baseMean,time2] = spikeZscoreSigTest( R2NoSound_Laser, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
            [R2NoSound_laser_Z, spikeMat, time1, indBase] = makeSpikeRates(R2NoSound_Laser, range , bin, filtWidth);
            
            subplot(3,4,[11,12]);
            plot( time1, R2NoSound_laser_Z,'color',laserRed,'linewidth',3);
            box off; set(gca,'tickdir','out','fontsize',16);
            xlim(xplotLim);
            yLims = get(gca, 'YLim');
            line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        end;
        
        %
        print(PSTH_Raster_plot, [FigSaveDir 'MD_Cell_' D '_Mixed_' num2str(i)], '-dpng');
        close all
    catch err
    end
end

%%
%%
for i = 1:numel(Spfc_mix)
    
    
    try
        PSTH_Raster_plot =  figure(i);
        PSTH_Raster_plot.Renderer='Painters';
        set(PSTH_Raster_plot, 'position', [63  39  1871  1318])
        set(PSTH_Raster_plot,'color','w');
        set(0,'DefaultAxesFontSize',16)
        set(0,'DefaultLineLineWidth',1)
        hold on;
        
        R1Sound = Spfc_mix(i).SpikeTimes_R1_sound( goodT(Spfc_mix(i).SpikeTimes_R1_sound));
        R2Sound = Spfc_mix(i).SpikeTimes_R2_sound( goodT(Spfc_mix(i).SpikeTimes_R2_sound));
        R1NoSound = Spfc_mix(i).SpikeTimes_R1_nosound( goodT(Spfc_mix(i).SpikeTimes_R1_nosound));
        R2NoSound = Spfc_mix(i).SpikeTimes_R2_nosound( goodT(Spfc_mix(i).SpikeTimes_R2_nosound));
        
        R1Sound_Laser = [ Spfc_mix(i).SpikeTimes_R1_sound_laser, Spfc_mix(i).SpikeTimes_R1_sound_inc_laser];
        R2Sound_Laser = [ Spfc_mix(i).SpikeTimes_R2_sound_laser, Spfc_mix(i).SpikeTimes_R2_sound_inc_laser];
        R1NoSound_Laser = [ Spfc_mix(i).SpikeTimes_R1_nosound_laser, Spfc_mix(i).SpikeTimes_R1_nosound_inc_laser];
        R2NoSound_Laser = [ Spfc_mix(i).SpikeTimes_R2_nosound_laser, Spfc_mix(i).SpikeTimes_R2_nosound_inc_laser];
        
        subplot(3,4,1);
        plotRaster( R1Sound );
        xlim(xplotLim);
        title('R1 - Aud','color','k');
        ylim( yplotLim );
        
        subplot(3,4,2);
        plotRaster( R2Sound );
        xlim(xplotLim);
        ylim( yplotLim );
        title('R2 - Aud','color','k');
        
        subplot(3,4,3);
        plotRaster( R1NoSound );
        xlim(xplotLim);
        ylim( yplotLim );
        title('R1 - Vis','color','r');
        
        subplot(3,4,4);
        plotRaster( R2NoSound );
        xlim(xplotLim);
        ylim( yplotLim );
        title('R2 - Vis','color','r');
        
        if size(R1Sound_Laser,2)~= 0
            subplot(3,4,5);
            plotRaster( R1Sound_Laser );
            xlim(xplotLim);
            title('R1 - Aud (Laser)','color',laserBlack);
            ylim( yplotLim );
        end;
        
        if size(R2Sound_Laser,2)~= 0
            subplot(3,4,6);
            plotRaster( R2Sound_Laser );
            xlim(xplotLim);
            ylim( yplotLim );
            title('R2 - Aud (Laser)','color',laserBlack);
        end;
        
        if size(R1NoSound_Laser,2)~= 0
            subplot(3,4,7);
            plotRaster( R1NoSound_Laser );
            xlim(xplotLim);
            ylim( yplotLim );
            title('R1 - Vis (Laser)','color',laserRed);
        end
        
        if size(R2NoSound_Laser)~=0
            subplot(3,4,8);
            plotRaster( R2NoSound_Laser );
            xlim(xplotLim);
            ylim( yplotLim );
            title('R2 - Vis (Laser)','color',laserRed);
        end
        
        % R1 - Aud ============================================================
        [R1SoundZ,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  R1Sound, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
        
        %         [R1SoundZ, spikeMat, time1, indBase] = makeSpikeRates(R1Sound, range , bin, filtWidth);
        
        subplot(3,4,[9,10]);
        plot( time1, R1SoundZ,'color','k','linewidth',3); hold on;
        box off; set(gca,'tickdir','out','fontsize',16);
        xlim(xplotLim);
        yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        title('Rule 1');
        
        % R1 - vis ============================================================
        [R1NoSoundZ,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  R1NoSound, 0.02, [-0.2 0.0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
        %         [R1NoSoundZ, spikeMat, time1, indBase] = makeSpikeRates(R1NoSound, range , bin, filtWidth);
        
        
        subplot(3,4,[9,10]);
        plot( time1, R1NoSoundZ,'color','r','linewidth',3);
        box off; set(gca,'tickdir','out','fontsize',16);
        xlim(xplotLim);     yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        
        if size(R1Sound_Laser,2)~= 0
            % R1 - Aud laser ============================================================
            %     [R1Sound_laser_Z,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  R1Sound_Laser, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
            [R1Sound_laser_Z, spikeMat, time1, indBase] = makeSpikeRates(R1Sound_Laser, range , bin, filtWidth);
            
            
            subplot(3,4,[9,10]);
            plot( time1, R1Sound_laser_Z,'color',laserBlack,'linewidth',3); hold on;
            box off; set(gca,'tickdir','out','fontsize',16);
            xlim(xplotLim);
            yLims = get(gca, 'YLim');
            line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        end
        
        if size(R1NoSound_Laser,2)~= 0
            % R1 - vis laser ============================================================
            %     [R1NoSound_laser_Z,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  R1NoSound_Laser, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
            [R1NoSound_laser_Z, spikeMat, time1, indBase] = makeSpikeRates(R1NoSound_Laser, range , bin, filtWidth);
            
            subplot(3,4,[9,10]);
            plot( time1, R1NoSound_laser_Z,'color',laserRed,'linewidth',3);
            box off; set(gca,'tickdir','out','fontsize',16);
            xlim(xplotLim);
            yLims = get(gca, 'YLim');
            line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        end;
        
        
        
        % R2 - Aud ============================================================
        %             [R2SoundZ,spikeMatFilt2,baseStd,baseMean,time2] = spikeZscoreSigTest( R2Sound, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
        [R2SoundZ, spikeMat, time1, indBase] = makeSpikeRates(R2Sound, range , bin, filtWidth);
        
        
        subplot(3,4,[11,12]);
        plot( time1, R2SoundZ,'color','k','linewidth',3);hold on;
        box off; set(gca,'tickdir','out','fontsize',16);
        xlim(xplotLim);
        yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        title('Rule 2');
        
        % R2 - Vis ============================================================
        %              [R2NoSoundZ,spikeMatFilt2,baseStd,baseMean,time2] = spikeZscoreSigTest( R2NoSound, 0.02,[-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
        [R2NoSoundZ, spikeMat, time1, indBase] = makeSpikeRates(R2NoSound, range , bin, filtWidth);
        
        
        subplot(3,4,[11,12]);
        plot( time1, R2NoSoundZ,'color','r','linewidth',3);
        box off; set(gca,'tickdir','out','fontsize',16);
        xlim(xplotLim);
        yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        
        
        if size(R2Sound_Laser,2)~= 0
            % R2 - Aud - laser ============================================================
            %     [R2Sound_laser_Z,spikeMatFilt2,baseStd,baseMean,time2] = spikeZscoreSigTest( R2Sound_Laser, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
            [R2Sound_laser_Z, spikeMat, time1, indBase] = makeSpikeRates(R2Sound_Laser, range , bin, filtWidth);
            
            subplot(3,4,[11,12]);
            plot( time1, R2Sound_laser_Z,'color',laserBlack,'linewidth',3);hold on;
            box off; set(gca,'tickdir','out','fontsize',16);
            xlim(xplotLim);
            yLims = get(gca, 'YLim');
            line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
            
        end;
        
        if size(R2NoSound_Laser,2)~= 0
            % R2 - Vis - laser ============================================================
            %     [R2NoSound_laser_Z,spikeMatFilt2,baseStd,baseMean,time2] = spikeZscoreSigTest( R2NoSound_Laser, 0.02, [-0.2 0], -0.5, 1.7 ,0.1,[], 0.5, 1.5 );
            [R2NoSound_laser_Z, spikeMat, time1, indBase] = makeSpikeRates(R2NoSound_Laser, range , bin, filtWidth);
            
            subplot(3,4,[11,12]);
            plot( time1, R2NoSound_laser_Z,'color',laserRed,'linewidth',3);
            box off; set(gca,'tickdir','out','fontsize',16);
            xlim(xplotLim);
            yLims = get(gca, 'YLim');
            line([0,0], [yLims]); line([preCueEnd,preCueEnd], [yLims]); line([delayEnd,delayEnd], [yLims]);
        end;
        
        %
        print(PSTH_Raster_plot, [FigSaveDir 'PFC_Cell_' D '_Mixed_' num2str(i)], '-dpng');
        close all
    catch err
    end
end

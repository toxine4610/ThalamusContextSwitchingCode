
clear all;
A = load('/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/DoubleCueDatabase/SOMCre/Learning/2018-03-21/SOMCre_Session.mat');
B = load('/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/DoubleCueDatabase/SOMCre/Learning/2018-03-21/SOMCre_Session_Change.mat');
%%
FigSaveDir = [D '_BPM33/'];
mkdir(FigSaveDir);
ANA = 0;

laserBlack = [0.5, 0.5, 0.5];
laserRed   = [219./255, 112./255, 147./255];

%%
indCorrVis = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 1);
indCorrAud = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 2);

Zaud = Z_C1( indCorrAud, :);
Zvis = Z_C1( indCorrVis, :);

indAudCS1 = find( Zaud(:,9) == 0);
indVisCS1 = find( Zvis(:,9) == 0);

indAudAssoc = find( Zaud(:,9) == 3);
indVisAssoc = find( Zvis(:,9) == 3);

indAudPT = find( Zaud(:,9) == 5);
indVisPT = find( Zvis(:,9) == 5);

[~, ~, goodPFC, goodMD] = cleanData(Spfc, Smd, Z_C1);
Spfc =  Spfc(goodPFC == 1);
Smd  =  Smd(goodMD == 1);


%%

delayEnd = 0.8 ; preCueEnd = 0.2;

range = [ -0.4, 1.5 ];
bin = 0.0080;
filtWidth = 0.1;

xplotLim = [-0.25, 1.2];
yplotLim = [0,20];

learn = 1;

%%

for i = 1:numel(Smd)
   
    PSTH_Raster_plot =  figure(i);
    PSTH_Raster_plot.Renderer='Painters';
    set(PSTH_Raster_plot, 'position', [63  39  1871  1318])
    set(PSTH_Raster_plot,'color','w');
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultLineLineWidth',1)
    hold on;
    
    Vis_CS1 = Smd(i).SpikeTimes_R1C1( indVisCS1);
    Aud_CS1 = Smd(i).SpikeTimes_R2C1( indAudCS1);
    
    Vis_CS1 = Vis_CS1( removeNoiseBursts(Vis_CS1) );
    Aud_CS1 = Aud_CS1( removeNoiseBursts( Aud_CS1) );
    
    Vis_Assoc = Smd(i).SpikeTimes_R1C1( indVisAssoc);
    Aud_Assoc = Smd(i).SpikeTimes_R2C1( indAudAssoc);
    
    Vis_Assoc = Vis_Assoc( removeNoiseBursts(Vis_Assoc) );
    Aud_Assoc = Aud_Assoc( removeNoiseBursts( Aud_Assoc) );
    
    if learn == 1
        Vis_CS2 = Smd(i).SpikeTimes_R1C1( indVisPT );
        Aud_CS2 = Smd(i).SpikeTimes_R2C1( indAudPT );
    end;
    
    subplot(3,3,1);
    plotRaster( Vis_CS1 );
    xlim(xplotLim); ylim( yplotLim );
    title('R1 (HP)','color','k');
    
    subplot(3,3,4);
    plotRaster( Aud_CS1 );
    xlim(xplotLim); ylim( yplotLim );
    title('R2 (LP)','color','r');
    
    subplot(3,3,2);
    plotRaster( Vis_Assoc );
    xlim(xplotLim); ylim( yplotLim );
    title('R1 (UV)','color',laserBlack);
    
    subplot(3,3,5);
    plotRaster( Aud_Assoc );
    xlim(xplotLim); ylim( yplotLim );
    title('R2 (Green)','color',laserRed);
    
    if learn == 1
        subplot(3,3,3);
        plotRaster( Vis_CS2 );
        xlim(xplotLim); ylim( yplotLim );
        title('R1 (Pure Tone)','color','k');
        
        subplot(3,3,6);
        plotRaster( Aud_CS2 );
        xlim(xplotLim); ylim( yplotLim );
        title('R2 (Pure Tone)','color','r');
    end;
    
    [Vis_CS1_rate, spikeMat, time1, indBase] = makeSpikeRates(Vis_CS1, range , bin, filtWidth);
    
    [Vis_Assoc_rate, spikeMat, time1, indBase] = makeSpikeRates(Vis_Assoc, range , bin, filtWidth);
    
    [Aud_CS1_rate, spikeMat, time1, indBase] = makeSpikeRates(Aud_CS1, range , bin, filtWidth);
    [Aud_Assoc_rate, spikeMat, time1, indBase] = makeSpikeRates(Aud_Assoc, range , bin, filtWidth);
    
    if learn == 1
        [Vis_CS2_rate, spikeMat, time1, indBase] = makeSpikeRates(Vis_CS2, range , bin, filtWidth);
        [Aud_CS2_rate, spikeMat, time1, indBase] = makeSpikeRates(Aud_CS2, range , bin, filtWidth);
    end;
    
    subplot(3,3,7);
    plot( time1, (Vis_CS1_rate-mean(Vis_CS1_rate(indBase)))./std(Vis_CS1_rate(indBase) ),'color','k');hold on
    plot( time1, (Aud_CS1_rate-mean(Aud_CS1_rate(indBase)))./std(Aud_CS1_rate(indBase) ),'color','r');
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
    xlim(xplotLim);
    yLims = get(gca, 'YLim');
    line([0,0], [yLims]); line([delayEnd,delayEnd], [yLims]);
    
    
    subplot(3,3,8);
    plot( time1, (Vis_Assoc_rate-mean(Vis_Assoc_rate(indBase)))./std(Vis_Assoc_rate(indBase) ),'color',laserBlack);hold on
    plot( time1, (Aud_Assoc_rate-mean(Aud_Assoc_rate(indBase)))./std(Aud_Assoc_rate(indBase) ),'color',laserRed);
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
    xlim(xplotLim);
    yLims = get(gca, 'YLim');
    line([0,0], [yLims]); line([delayEnd,delayEnd], [yLims]);
    
    if learn == 1
        subplot(3,3,9);
        plot( time1, (Vis_CS2_rate-mean(Vis_CS2_rate(indBase)))./std(Vis_CS2_rate(indBase) ),'color','k');hold on
        plot( time1, (Aud_CS2_rate-mean(Aud_CS2_rate(indBase)))./std(Aud_CS2_rate(indBase) ),'color','r');
        axis normal; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
        xlim(xplotLim);
        yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([delayEnd,delayEnd], [yLims]);
    end
%     pause
    print(PSTH_Raster_plot, [FigSaveDir 'MD_Cell_' D '_Learning_' num2str(i)], '-dpng');
    close all
end;

%%

for i = 1:numel(Spfc)
    
    PSTH_Raster_plot =  figure(i);
    PSTH_Raster_plot.Renderer='Painters';
    set(PSTH_Raster_plot, 'position', [63  39  1871  1318])
    set(PSTH_Raster_plot,'color','w');
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultLineLineWidth',1)
    hold on;
    
    Vis_CS1 = Smd(i).SpikeTimes_R1C1(indVisCS1);
    Aud_CS1 = Smd(i).SpikeTimes_R2C1(indAudCS1);
    
    Vis_CS1 = Vis_CS1( removeNoiseBursts(Vis_CS1) );
    Aud_CS1 = Aud_CS1( removeNoiseBursts( Aud_CS1) );
    
    Vis_Assoc = Smd(i).SpikeTimes_R1C1( indVisAssoc);
    Aud_Assoc = Smd(i).SpikeTimes_R2C1( indAudAssoc);
    
    Vis_Assoc = Vis_Assoc( removeNoiseBursts(Vis_Assoc) );
    Aud_Assoc = Aud_Assoc( removeNoiseBursts( Aud_Assoc) );
    
    
    if learn == 1
        Vis_CS2 = Spfc(i).SpikeTimes_R1C1( indVisPT ) ;
        Aud_CS2 = Spfc(i).SpikeTimes_R2C1( indAudPT );
    end;
    
    subplot(3,3,1);
    plotRaster( Vis_CS1 );
    xlim(xplotLim); ylim( yplotLim );
    title('R1 (HP)','color','k');
    
    
    subplot(3,3,4);
    plotRaster( Aud_CS1 );
    xlim(xplotLim); ylim( yplotLim );
    title('R2 (LP)','color','r');
    
    subplot(3,3,2);
    plotRaster( Vis_Assoc );
    xlim(xplotLim); ylim( yplotLim );
    title('R1 (UV)','color',laserBlack);
    
    subplot(3,3,5);
    plotRaster( Aud_Assoc );
    xlim(xplotLim); ylim( yplotLim );
    title('R2 (Green)','color',laserRed);
    
    if learn == 1
        subplot(3,3,3);
        plotRaster( Vis_CS2 );
        xlim(xplotLim); ylim( yplotLim );
        title('R1 (1pip @9Khz)','color','k');
        
        subplot(3,3,6);
        plotRaster( Aud_CS2 );
        xlim(xplotLim); ylim( yplotLim );
        title('R2 (3pips @9Khz)','color','r');
    end;
    
    
    [Vis_CS1_rate, spikeMat, time1, indBase] = makeSpikeRates(Vis_CS1, range , bin, filtWidth);       
%     [Vis_CS2_rate, spikeMat, time1, indBase] = makeSpikeRates(Vis_CS2, range , bin, filtWidth);
    [Vis_Assoc_rate, spikeMat, time1, indBase] = makeSpikeRates(Vis_Assoc, range , bin, filtWidth);
    
%     [Aud_CS1_rate, spikeMat, time1, indBase] = makeSpikeRates(Aud_CS1, range , bin, filtWidth);
    %         [Aud_CS2_rate, spikeMat, time1, indBase] = makeSpikeRates(Aud_CS2, range , bin, filtWidth);
    [Aud_Assoc_rate, spikeMat, time1, indBase] = makeSpikeRates(Aud_Assoc, range , bin, filtWidth);
    
    subplot(3,3,7);
    plot( time1, (Vis_CS1_rate-mean(Vis_CS1_rate(indBase)))./std(Vis_CS1_rate(indBase) ),'color','k');hold on
    plot( time1, (Aud_CS1_rate-mean(Aud_CS1_rate(indBase)))./std(Aud_CS1_rate(indBase) ),'color','r');
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
    xlim(xplotLim);
    yLims = get(gca, 'YLim');
    line([0,0], [yLims]); line([delayEnd,delayEnd], [yLims]);
    
    
    subplot(3,3,8);
    plot( time1, (Vis_Assoc_rate-mean(Vis_Assoc_rate(indBase)))./std(Vis_Assoc_rate(indBase) ),'color',laserBlack);hold on
    plot( time1, (Aud_Assoc_rate-mean(Aud_Assoc_rate(indBase)))./std(Aud_Assoc_rate(indBase) ),'color',laserRed);
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
    xlim(xplotLim);
    yLims = get(gca, 'YLim');
    line([0,0], [yLims]); line([delayEnd,delayEnd], [yLims]);
    
    
    if learn == 1
        subplot(3,3,9);
        plot( time1, (Vis_CS2_rate-mean(Vis_CS2_rate(indBase)))./std(Vis_CS2_rate(indBase) ),'color','k');hold on
        plot( time1, (Aud_CS2_rate-mean(Aud_CS2_rate(indBase)))./std(Aud_CS2_rate(indBase) ),'color','r');
        axis normal; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
        xlim(xplotLim);
        yLims = get(gca, 'YLim');
        line([0,0], [yLims]); line([delayEnd,delayEnd], [yLims]);
    end
    
    print(PSTH_Raster_plot, [FigSaveDir 'PFC_Cell_' D '_Learning_' num2str(i)], '-dpng');
    close all
end;
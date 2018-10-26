function Shrew_plotFigures(D, Z_C1, Spfc, Smd)

hasMix = 0;

FigSaveDir = ['../Shrew_Figrs/' D '/'];

% FigSaveDir = [ 'C:\Users\Halassalab-CG\Dropbox\Rajeev\Shrew Peach\2018-01-26_11-37-16\' D '\'];
mkdir(FigSaveDir);

binS = 0.005;FW = 0.1;

% indAud = find( Z_C1(:,4)==80 & abs(Z_C1(:,5))==2 & Z_C1(:,1)==1 );
% indVis = find( Z_C1(:,4)==80 & abs(Z_C1(:,5))==1 & Z_C1(:,1)==1 );
% 
% Zaud  =  Z_C1(indAud, :);
% Zvis  =  Z_C1(indVis, :);
% 
% indAud = find(Zaud(:,2) >= 0.60);
% indVis = find(Zvis(:,2) >= 0.60);

range = [ -0.2, 1.6 ];

bin = 0.0010;
filtWidth = 0.12;

%%

for i = 1:numel(Spfc)
    
    PSTH_Raster_plot =  figure(i);
    PSTH_Raster_plot.Renderer='Painters';
    set(PSTH_Raster_plot, 'position', [ 1         120        1280         636])
    set(PSTH_Raster_plot,'color','w');
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultLineLineWidth',1)
    hold on;
    
    figure(i); subplot(3,2,1);
    plotRaster(Spfc(i).SpikeTimes_Pure_Audition);xlim([-0.25,1.8]);
    title('Audition','color','r');
    
    figure(i); subplot(3,2,3);
    plotRaster(Spfc(i).SpikeTimes_Pure_Vision);xlim([-0.25,1.8]);
    title('Vision','color', 'b');
    
    if hasMix == 1
        figure(i); subplot(3,2,2);
        plotRaster(Spfc(i).SpikeTimes_Mix_Audition);xlim([-0.25,1.8]);
        title('Audition (60/40) Dropout','color','m');
        
        figure(i); subplot(3,2,4);
        plotRaster(Spfc(i).SpikeTimes_Mix_Vision);xlim([-0.25,1.8]);
        title('Vision (60/40) Dropout','color','c');
        
        [AM, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Mix_Audition, range , bin, filtWidth);
        [VM, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Mix_Vision, range , bin, filtWidth);
        
    end
    
    [AP, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Pure_Audition, range , bin, filtWidth);
    [VP, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Pure_Vision, range , bin, filtWidth);
    
    
    %     [VP,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( Spfc(i).SpikeTimes_Pure_Vision(indVis) , binS, [-0.1 0.2], -0.5, 1.8 ,FW,[], 0.5, 1.5 );
    %     [AP,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( Spfc(i).SpikeTimes_Pure_Audition(indAud)  , binS, [-0.1 0.2], -0.5, 1.8 ,FW,[], 0.5, 1.5 );
    
    %     [VM,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(  ,binS, [-0.1 0.2], -0.5, 1.8 ,FW,[], 0.5, 1.5 );
    %     [AM,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest(   ,binS, [-0.1 0.2], -0.5, 1.8 ,FW,[], 0.5, 1.5 );
    
    figure(i); subplot(3,2,5);
    plot( time1, VP-nanmean(VP(indBase)),'color','b','linewidth',3);
    hold on;
    plot( time1, AP-nanmean(AP(indBase)),'color','r','linewidth',3); hold on;
    xlim([-0.25,1.8]);
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);
    yLims = get(gca, 'YLim');
    line([0.0,0.0], [yLims]);
    line([1,1], [yLims]);
%     line([0.725,0.725], [yLims]);
%     line([0.925,0.925], [yLims]);
%     line([1.5,1.5], [yLims]);
    
    if hasMix == 1
        figure(i); subplot(3,2,6);
        plot( time1, VM-nanmean(VM(indBase)),'color','c','linewidth',3);
        hold on;
        plot( time1, AM-nanmean(VM(indBase)),'color','m','linewidth',3); hold on;
        xlim([-0.25,1.8]);
        axis normal; box off; set(gca,'tickdir','out','fontsize',16);
        yLims = get(gca, 'YLim');
        line([0.0,0.0], [yLims]);
        line([0.5,0.5], [yLims]);
%         line([0.725,0.725], [yLims]);
%         line([0.925,0.925], [yLims]);
%         line([1.5,1.5], [yLims]);
    end;
    
    print(PSTH_Raster_plot, [FigSaveDir 'PFC_Cell_' D '_Accum_' num2str(i)], '-dpng');
    close all;
end;


%%

for i = 1:numel(Smd)
    
    PSTH_Raster_plot =  figure(i);
    PSTH_Raster_plot.Renderer='Painters';
    set(PSTH_Raster_plot, 'position', [ 1         120        1280         636])
    set(PSTH_Raster_plot,'color','w');
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultLineLineWidth',1)
    hold on;
    
    figure(i); subplot(3,2,1);
    plotRaster(Smd(i).SpikeTimes_Pure_Audition);xlim([-0.25,1.8]);
    title('Audition','color','r');
    
    figure(i); subplot(3,2,3);
    plotRaster(Smd(i).SpikeTimes_Pure_Vision);xlim([-0.25,1.8]);
    title('Vision ','color', 'b');
    
    if hasMix == 1
        figure(i); subplot(3,2,2);
        plotRaster(Smd(i).SpikeTimes_Mix_Audition);xlim([-0.25,1.8]);
        title('Audition (60/40) Dropout','color','m');
        
        figure(i); subplot(3,2,4);
        plotRaster(Smd(i).SpikeTimes_Mix_Vision);xlim([-0.25,1.8]);
        title('Vision (60/40) Dropout','color','c');
        
        [AM, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Mix_Audition, range , bin, filtWidth);
        [VM, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Mix_Vision, range , bin, filtWidth);
        
    end
    
    [AP, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Pure_Audition, range , bin, filtWidth);
    [VP, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Pure_Vision, range , bin, filtWidth);
    
    figure(i); subplot(3,2,5);
    plot( time1, VP,'color','b','linewidth',3);
    hold on;
    plot( time1, AP,'color','r','linewidth',3); hold on;
    xlim([-0.25,1.8]);
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);
    yLims = get(gca, 'YLim');
    line([0.0,0.0], [yLims]);
%     line([0.5,0.5], [yLims]);
%     line([0.725,0.725], [yLims]);
%     line([0.925,0.925], [yLims]);
    line([1.0,1.0], [yLims]);
    
    if hasMix == 1
        figure(i); subplot(3,2,6);
        plot( time1, VM,'color','c','linewidth',3);
        hold on;
        plot( time1, AM,'color','m','linewidth',3); hold on;
        xlim([-0.25,1.8]);
        axis normal; box off; set(gca,'tickdir','out','fontsize',16);
        yLims = get(gca, 'YLim');
        line([0.2,0.2], [yLims]);
        line([0.5,0.5], [yLims]);
        line([0.725,0.725], [yLims]);
        line([0.925,0.925], [yLims]);
        line([1.5,1.5], [yLims]);
    end;
    
    print(PSTH_Raster_plot, [FigSaveDir 'MD_Cell_' D '_Accum_' num2str(i)], '-dpng');
    close all;
end;
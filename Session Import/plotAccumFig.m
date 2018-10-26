FigSaveDir = ['../Acc1_Figrs/' D '/'];
mkdir(FigSaveDir);

%%

binS = 0.005;FW = 0.1;

indAud = find( Z_C1(:,4)==80 & abs(Z_C1(:,5))==2 & Z_C1(:,1)==1 );
indVis = find( Z_C1(:,4)==80 & abs(Z_C1(:,5))==1 & Z_C1(:,1)==1 );

indAud2 = find( Z_C1(:,4)==60 & abs(Z_C1(:,5))==2 & Z_C1(:,1)==1 );
indVis2 = find( Z_C1(:,4)==60 & abs(Z_C1(:,5))==1 & Z_C1(:,1)==1 );


Zaud  =  Z_C1(indAud, :);
Zvis  =  Z_C1(indVis, :);

Zaud2  =  Z_C1(indAud2, :);
Zvis2  =  Z_C1(indVis2, :);

indAud = find(Zaud(:,2) >= 0.50);
indVis = find(Zvis(:,2) >= 0.50);

indAudLaser = find(Zaud2(:,9) == 1 & Zaud2(:,4) == 60);
indVisLaser = find(Zvis2(:,9) == 1 & Zvis2(:,4) == 60);

indAudNoLaser = find(Zaud2(:,9) == 0 & Zaud2(:,4) == 60);
indVisNoLaser = find(Zvis2(:,9) == 0 & Zvis2(:,4) == 60);


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
    plotRaster(Spfc(i).SpikeTimes_Pure_Audition(indAud));xlim([-0.25,1.8]);
    title('Audition (60/40)','color','r');
    
    figure(i); subplot(3,2,3);
    plotRaster(Spfc(i).SpikeTimes_Pure_Vision(indVis));xlim([-0.25,1.8]);
    title('Vision (60/40)','color', 'b');
    
    figure(i); subplot(3,2,2);
    plotRaster(Spfc(i).SpikeTimes_Mix_Audition(indAudLaser));xlim([-0.25,1.8]);
    title('Audition Pure Laser','color','m');
    
    figure(i); subplot(3,2,4);
    plotRaster(Spfc(i).SpikeTimes_Mix_Vision(indVisLaser));xlim([-0.25,1.8]);
    title('Vision Pure Laser','color','c');
    
    [AM, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Mix_Audition(indAudLaser), range , bin, filtWidth);
    [VM, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Mix_Vision(indVisLaser), range , bin, filtWidth);
    
     
    [AMNL, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Mix_Audition(indAudNoLaser), range , bin, filtWidth);
    [VMNL, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Mix_Vision(indVisNoLaser), range , bin, filtWidth);
    
    [AP, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Pure_Audition(indAud), range , bin, filtWidth);
    [VP, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Spfc(i).SpikeTimes_Pure_Vision(indVis), range , bin, filtWidth);
    
    
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
    
    line([1.5,1.5], [yLims]);
    
    
    figure(i); subplot(3,2,6);
    plot( time1, VM-nanmean(VM(indBase)),'color','c','linewidth',3);
    hold on;
    plot( time1, AM-nanmean(VM(indBase)),'color','m','linewidth',3); hold on;
    
    plot( time1, VMNL-nanmean(AMNL(indBase)),'color','b','linewidth',3);
    hold on;
    plot( time1, AMNL-nanmean(AMNL(indBase)),'color','r','linewidth',3); hold on;
    
    xlim([-0.25,1.8]);
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);
    yLims = get(gca, 'YLim');
   
    line([1.5,1.5], [yLims]);
    
    
     durPulse = 50; durISI = 25;
    Npulse = 9;
    tvec = [0:(durPulse+durISI):Npulse*(durPulse+durISI)]/1000;
    for j = 1:length(tvec);
        figure(i); subplot(3,2,6);
        yLims = get(gca, 'YLim');
        line([tvec(j),tvec(j)], [yLims]);
        
        figure(i); subplot(3,2,5);
        yLims = get(gca, 'YLim');
        line([tvec(j),tvec(j)], [yLims]);
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
    plotRaster(Smd(i).SpikeTimes_Pure_Audition(indAud));xlim([-0.25,1.8]);
    title('Audition (60/40)','color','r');
    ylim([0,20]);
    
    figure(i); subplot(3,2,3);
    plotRaster(Smd(i).SpikeTimes_Pure_Vision(indVis));xlim([-0.25,1.8]);
    title('Vision (60/40)','color', 'b');
     ylim([0,20]);
     
    figure(i); subplot(3,2,2);
    plotRaster(Smd(i).SpikeTimes_Mix_Audition(indAudLaser));xlim([-0.25,1.8]);
    title('Audition Laser','color','m');
     ylim([0,20]);
     
    figure(i); subplot(3,2,4);
    plotRaster(Smd(i).SpikeTimes_Mix_Vision(indVisLaser));xlim([-0.25,1.8]);
    title('Vision Laser','color','c');
     ylim([0,20]);
     
    [AMNL, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Mix_Audition(indAudNoLaser), range , bin, filtWidth);
    [VMNL, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Mix_Vision(indVisNoLaser), range , bin, filtWidth);
    
     
    [AM, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Mix_Audition(indAudLaser), range , bin, filtWidth);
    [VM, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Mix_Vision(indVisLaser), range , bin, filtWidth);
    
    [AP, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Pure_Audition(indAud), range , bin, filtWidth);
    [VP, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(Smd(i).SpikeTimes_Pure_Vision(indVis), range , bin, filtWidth);
    
    figure(i); subplot(3,2,5);
    plot( time1, VP,'color','b','linewidth',3);
    hold on;
    plot( time1, AP,'color','r','linewidth',3); hold on;
    xlim([-0.25,1.8]);  ylim([0,2]);
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);
    yLims = get(gca, 'YLim');
    line([0,0.0], [yLims]);
    line([0.725,0.725], [yLims]);
    line([1.5,1.5], [yLims]);
    
    
    figure(i); subplot(3,2,6);
    plot( time1, VM,'color','c','linewidth',3);
    hold on;
    plot( time1, AM,'color','m','linewidth',3); hold on;
    
    plot( time1, VMNL,'color','b','linewidth',3);
    hold on;
    plot( time1, AMNL,'color','r','linewidth',3); hold on;
    
    xlim([-0.25,1.8]);  ylim([0,2]);
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);
    yLims = get(gca, 'YLim');
    line([0,0.0], [yLims]);
    line([0.725,0.725], [yLims]);
    line([1.5,1.5], [yLims]);
    
    durPulse = 50; durISI = 25;
    Npulse = 9;
    tvec = [0:(durPulse+durISI):Npulse*(durPulse+durISI)]/1000;
    for j = 1:length(tvec);
        figure(i); subplot(3,2,6);
        yLims = get(gca, 'YLim');
        line([tvec(j),tvec(j)], [yLims]);
        
        figure(i); subplot(3,2,5);
        yLims = get(gca, 'YLim');
        line([tvec(j),tvec(j)], [yLims]);
    end;


        
    
    print(PSTH_Raster_plot, [FigSaveDir 'MD_Cell_' D '_Accum_' num2str(i)], '-dpng');
    close all;
end;
FigSaveDir = ['C:\Users\Halassalab-CG\Dropbox\Rajeev\ForContextSwitchProjectt\Mixed_SOMCRE_Figrs\' D '\'];
mkdir(FigSaveDir);
%%

for i = 1:numel(P1.Smd_Mixt)
    
    R1Sound = [ P1.Smd_Mixt(i).R1Sound, P2.Smd_Mixt(i).R1Sound ];
    R2Sound = [ P1.Smd_Mixt(i).R2Sound, P2.Smd_Mixt(i).R2Sound ];
    R1NoSound = [ P1.Smd_Mixt(i).R1NoSound, P2.Smd_Mixt(i).R1NoSound ];
    R2NoSound = [ P1.Smd_Mixt(i).R2NoSound, P2.Smd_Mixt(i).R2NoSound ];
    
    
    PSTH_Raster_plot =  figure(i);
    set(PSTH_Raster_plot, 'position', [43         188        1041         517])
    set(PSTH_Raster_plot,'color','w');
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultLineLineWidth',1)
    hold on;
    
    figure(i); set(gcf,'color','w')
    figure(i);subplot(3,2,1);
    plotRaster( R1Sound ); xlim([-.3, 1.2]);
    title('R1 Sound');
    
    figure(i);subplot(3,2,2);
    plotRaster( R2Sound ); xlim([-.3, 1.2]);
    title('R2 Sound');
    
    
    figure(i);subplot(3,2,3);
    plotRaster( R1NoSound ); xlim([-.3, 1.2]);
    title('R1 No Sound');
    
    figure(i);subplot(3,2,4);
    plotRaster( R2NoSound ); xlim([-.3, 1.2]);
    title('R2 No Sound');
    
    figure(i);subplot( 3,2,[5]);
    
    [sound,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( R1Sound , 0.025, [0.0 0.25], -0.5, 1.2 ,0.1,[], 0.5, 1.5 );
    [nosound,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( R1NoSound , 0.025, [0.0 0.25], -0.5, 1.2 ,0.1,[], 0.5, 1.5 );
    plot( time1, sound,'k','linewidth',3); hold on;
    plot( time1, nosound, 'r','linewidth',3);
    xlim([-.3, 1.2]);
    yLims = get(gca, 'YLim');
    line([0,0], [yLims]);
    line([0.3,0.3], [yLims]);
    line([0.9,0.9], [yLims]);
    box off;
    
    figure(i);subplot( 3,2,[6]);
    
    [sound,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( R2Sound , 0.025, [0.0 0.25], -0.5, 1.2 ,0.1,[], 0.5, 1.5 );
    [nosound,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( R2NoSound , 0.025, [0.0 0.25], -0.5, 1.2 ,0.1,[], 0.5, 1.5 );
    plot( time1, sound,'k','linewidth',3); hold on;
    plot( time1, nosound, 'r','linewidth',3);
    xlim([-.3, 1.2]);
    yLims = get(gca, 'YLim');
    line([0,0], [yLims]);
    line([0.3,0.3], [yLims]);
    line([0.9,0.9], [yLims]);
    box off;
    
    print(PSTH_Raster_plot, [FigSaveDir 'MD_Cell_MixCue_' num2str(i)], '-dpng');
    close all;
    
end;

%%

for i = 47
    
    R1Sound = [ P1.Spfc_Mixt(i).R1Sound, P2.Spfc_Mixt(i).R1Sound ];
    R2Sound = [ P1.Spfc_Mixt(i).R2Sound, P2.Spfc_Mixt(i).R2Sound ];
    R1NoSound = [ P1.Spfc_Mixt(i).R1NoSound, P2.Spfc_Mixt(i).R1NoSound ];
    R2NoSound = [ P1.Spfc_Mixt(i).R2NoSound, P2.Spfc_Mixt(i).R2NoSound ];
    
    
    PSTH_Raster_plot =  figure(i);
    set(PSTH_Raster_plot, 'position', [43         188        1041         517])
    set(PSTH_Raster_plot,'color','w');
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultLineLineWidth',1)
    hold on;
    
    figure(i); set(gcf,'color','w')
    figure(i);subplot(3,2,1);
    plotRaster( R1Sound ); xlim([-.3, 1.2]);
    title('R1 Sound');
    
    figure(i);subplot(3,2,2);
    plotRaster( R2Sound ); xlim([-.3, 1.2]);
    title('R2 Sound');
    
    
    figure(i);subplot(3,2,3);
    plotRaster( R1NoSound ); xlim([-.3, 1.2]);
    title('R1 No Sound');
    
    figure(i);subplot(3,2,4);
    plotRaster( R2NoSound ); xlim([-.3, 1.2]);
    title('R2 No Sound');
    
    figure(i);subplot( 3,2,[5]);
    
    [sound,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( R1Sound , 0.02, [0.0 0.3], -0.5, 1.2 ,0.1,[], 0.5, 1.5 );
    [nosound,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( R1NoSound , 0.02, [0.0 0.3], -0.5, 1.2 ,0.1,[], 0.5, 1.5 );
    plot( time1, sound,'k','linewidth',3); hold on;
    plot( time1, nosound, 'r','linewidth',3);
    xlim([-.3, 1.2]);
    yLims = get(gca, 'YLim');
    line([0,0], [yLims]);
    line([0.3,0.3], [yLims]);
    line([0.9,0.9], [yLims]);
    box off;
    
    figure(i);subplot( 3,2,[6]);
    
    [sound,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( R2Sound , 0.02, [0.0 0.3], -0.5, 1.2 ,0.1,[], 0.5, 1.5 );
    [nosound,spikeMatFilt1,baseStd,baseMean,time1] = spikeZscoreSigTest( R2NoSound , 0.02, [0.0 0.3], -0.5, 1.2 ,0.1,[], 0.5, 1.5 );
    plot( time1, sound,'k','linewidth',3); hold on;
    plot( time1, nosound, 'r','linewidth',3);
    xlim([-.3, 1.2]);
    yLims = get(gca, 'YLim');
    line([0,0], [yLims]);
    line([0.3,0.3], [yLims]);
    line([0.9,0.9], [yLims]);
    box off;
    
%     print(PSTH_Raster_plot, [FigSaveDir 'PFC_Cell_MixCue_' num2str(i)], '-dpng');
%     close all;
    
end;
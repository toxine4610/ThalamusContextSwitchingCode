PFCTT = 1:16;
MDTT  = 17:32;

ct1 = 0;
ct2 = 0;

range = [ -0.2, 1.6 ];

bin = 0.0010;
filtWidth = 0.12;

for n = 1:size(num_seq,1)
    tetrode_num = num_seq( n,1 );
    
    if ismember( tetrode_num ,PFCTT)
        ct1 = ct1+1;
        
        Spfc(ct1).unitNbr = n;
        Spfc(ct1).ttNbr   = tetrode_num;
        eval( ['foo = cl' num2str(n) '_cellDatabaseline;'] );
        if isfield(foo, 'Vis_corr')
            SpikeTimes_VC = foo.Vis_corr;
        else
            SpikeTimes_VC = [];
        end;
        
        if isfield(foo, 'Aud_corr')
            SpikeTimes_AC = foo.Aud_corr;
        else
            SpikeTimes_AC = [];
        end;
        
        Spfc(ct1).SpikeTimes_Pure_Audition = [SpikeTimes_AC];
        Spfc(ct1).SpikeTimes_Pure_Vision   = [SpikeTimes_VC];
        clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
        
    elseif ismember( tetrode_num, MDTT)
        
        ct2 = ct2+1;
        
        Smd(ct2).unitNbr = n;
        Smd(ct2).ttNbr   = tetrode_num;
        
        eval( ['foo = cl' num2str(n) '_cellDatabaseline;'] );
        clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
        
        if isfield(foo, 'Vis_corr')
            SpikeTimes_VC = foo.Vis_corr;
        else
            SpikeTimes_VC = [];
        end;
        
        
        
        if isfield(foo, 'Aud_corr')
            SpikeTimes_AC = foo.Aud_corr;
        else
            SpikeTimes_AC = [];
        end;
        
        Smd(ct2).SpikeTimes_Pure_Audition = [SpikeTimes_AC];
        Smd(ct2).SpikeTimes_Pure_Vision   = [SpikeTimes_VC];
        clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
        
    end;
    
end;

%%







figure(1); subplot(3,2,1);
plotRaster( cl3_cellDatabaseline.Aud_corr );xlim([-0.25,1.8]);
title('Audition','color','r');

figure(1); subplot(3,2,3);
plotRaster( cl3_cellDatabaseline.Vis_corr );xlim([-0.25,1.8]);
title('Vision','color', 'b');


[AP, spikeMat_aud_mix, time1, indBase] = makeSpikeRatesSI(cl3_cellDatabaseline.Aud_corr , range , bin, filtWidth);
[VP, spikeMat_vis_mix, time1, indBase] = makeSpikeRatesSI(cl3_cellDatabaseline.Vis_corr , range , bin, filtWidth);

figure(1); subplot(3,2,5);

plot( time1, VP-nanmean(VP(indBase)),'color','b','linewidth',3);
hold on;
plot( time1, AP-nanmean(AP(indBase)),'color','r','linewidth',3); hold on;
xlim([-0.25,1.8]);
axis normal; box off; set(gca,'tickdir','out','fontsize',16);
yLims = get(gca, 'YLim');
line([1.0,1.0], [yLims]);
line([0,0], [yLims]);
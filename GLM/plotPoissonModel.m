
ind = find( time1 >= -0.2 & time1 <= 1.0 );
Rate = R1NoSoundZ( ind );
% Rate = Rate./max(Rate);

timeaxis = time1( ind );

dt = 0.5; % 0.02;

S = [];

for trial = 1:numel(R1NoSound)
    spikes = [];
    for i = 1:length(Rate)
        spikes(1, i) = (Rate(i).*dt ) > rand(1,1);
    end;
    S = cat(1, S, spikes);
end;

ModelRate =  smooth( mean(S,1), 20  );

figure(1);
yyaxis right;
plot( time1(ind), ModelRate./max(ModelRate) ); hold on

yyaxis left;
plot( time1(ind), Rate );
% 

 PSTH_Raster_plot =  figure(2);
 PSTH_Raster_plot.Renderer='Painters';

subplot(1,2,2);
plotRaster( R1NoSound );
xlim([-0.2, 1]);

subplot(1,2,1);
for trial = 1:20
    for t = 1:length(timeaxis)
        if S(trial,t)==0
            line( [timeaxis(t), timeaxis(t)], [trial-1, trial-1] ); hold on;
        elseif S(trial,t) == 1
            line( [timeaxis(t), timeaxis(t)], [trial-1, trial] ); hold on;
        end;
    end;
end;
xlim([-0.2, 1]); ylim([0,44]);

CMI   = @(X,Y) (X-Y)./(X+Y);
goodT =  @(X) find(~cellfun(@isempty,X));

Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/SOMCre/DualModalityCue/';
ExptDates = {'2017-12-15','2017-12-16','2017-12-17','2017-12-18','2017-12-20','2017-12-22','2017-12-23','2017-12-24','2017-12-26'};

range = [ -0.2, 1.6 ];
bin = 0.0010;
filtWidth = 0.1;
%%
for d = 1:length(ExptDates)
    
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} '/Somcre_mixed.mat'] )
    
    %%
    
    [Spfc_mix, Smd_mix, ZVis, ZAud] = packageData(Z_C1, Smd, Spfc);
    PVC1 = []; PVC2 = [];
    
    for i = 1:numel(Smd_mix)
        
        clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateR1C1 spikeRateR2C1 spikeRateR1C2 spikeRateR2C2
        % correct =============================================================
        R1C1  = Smd_mix(i).SpikeTimes_R1_sound;
        R2C1  = Smd_mix(i).SpikeTimes_R2_sound;
        
        R1C2  = Smd_mix(i).SpikeTimes_R1_nosound;
        R2C2  = Smd_mix(i).SpikeTimes_R2_nosound;
        
        R1C1  =  R1C1( goodT(R1C1) );
        R1C2  =  R1C2( goodT(R1C2) );
        R2C1  =  R2C1( goodT(R2C1) );
        R2C2  =  R2C2( goodT(R2C2) );
        
        
        C1 = [R1C1];
        C2 = [R2C1];
        
        [spikeRateC1, spikeMat_C1_mix, timeOut, indBase] = makeSpikeRates(C1, range , bin, filtWidth);
        [spikeRateC2, spikeMat_C2_mix, timeOut, indBase] = makeSpikeRates(C2, range , bin, filtWidth);
        
        R_C1(i,:) = ( spikeRateC1  - nanmin( spikeRateC1 ) ) ./ ( nanmax( spikeRateC1) - nanmin( spikeRateC1 ) );
        R_C2(i,:) = ( spikeRateC2  - nanmin( spikeRateC2 ) ) ./ ( nanmax( spikeRateC2) - nanmin( spikeRateC2 ) );
        
        indDelay = find( timeOut >=0.35 & timeOut <=0.85 );
        cmi(i)    = CMI( nanmean( R_C1(i, indDelay) ), nanmean( R_C2(i, indDelay)) );
        
        PVC1{i} = ( spikeMat_C1_mix );
        PVC2{i} = ( spikeMat_C2_mix );
        
    end;
    
    %%
    
    timeAccum = find( timeOut >= 0.35 & timeOut <= .85);
    
    
    [coeff,score,latent,tsquared,explained,mu] = pca(R_C1( find(cmi>0), timeAccum));
    
    [coeff_m,score_m,latent_m,tsquared_m,explained_m,mu_m] = pca(R_C2( find(cmi<0),timeAccum));
    
    figure(200); set(gcf,'color','w');
    subplot(1,2,1);
    plot( timeOut(timeAccum), coeff_m(:,1) ,'k'); hold on
    plot( timeOut(timeAccum), coeff(:,1),'r' ); hold on
    
    subplot(1,2,2);
    plot( timeOut(timeAccum), coeff_m(:,2) ,'k'); hold on
    plot( timeOut(timeAccum), coeff(:,2),'r' ); hold on
    
    figure(800);set(gcf,'color','w');
    plot3( smooth( coeff(:,1) ),smooth( coeff(:,2) ), smooth( coeff(:,3) ),'r'); hold on;
    plot3( smooth( coeff_m(:,1) ),smooth( coeff_m(:,2) ), smooth( coeff_m(:,3) ),'k'); hold on;
    
    plot3( coeff(1,1) ,smooth( coeff(1,2) ), smooth( coeff(1,3) ),'ko'); hold on;
    plot3( coeff(end,1) ,smooth( coeff(end,2) ), smooth( coeff(end,3) ),'sk'); hold on;
    
    plot3( coeff_m(1,1) ,smooth( coeff_m(1,2) ), smooth( coeff_m(1,3) ),'ro'); hold on;
    plot3( coeff_m(end,1) ,smooth( coeff_m(end,2) ), smooth( coeff_m(end,3) ),'sr'); hold on;
    
    Diff(d,:) = coeff(:,1) - coeff_m(:,1);
    
    %%
    timeDelay = find( timeOut >= 0.0 & timeOut <= .85);
    
    clear PVC1n PVC2n
    
    NTrials = cellfun( @(x) size(x,1), PVC1 );
    NTrialskeep = min(NTrials);
    
    for i = 1:numel(Smd_mix)
        this_cell = i;
        foo =  PVC1{this_cell}( 1:NTrialskeep, timeDelay );
        PVC1n(i,:,:) = foo;
    end;
    
    NTrials = cellfun( @(x) size(x,1), PVC2 );
    NTrialskeep = min(NTrials);
    
    for i = 1:numel(Smd_mix)
        this_cell = i;
        foo =  PVC2{this_cell}( 1:NTrialskeep, timeDelay );
        PVC2n(i,:,:) = foo;
    end;
    
    PVC1n = permute( PVC1n, [2,1,3] );
    PVC2n = permute( PVC2n, [2,1,3] );
    
    PVC1n = squeeze( nanmean( PVC1n, 3) );
    PVC2n = squeeze( nanmean( PVC2n, 3) );

    
end


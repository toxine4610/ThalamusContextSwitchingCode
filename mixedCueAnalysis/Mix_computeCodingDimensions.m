
%%
CD = @(r1, r2) ( nanmean(r1,1) - nanmean(r2,1) )./sqrt( nanvar(r1,[],1)+nanvar(r2,[],1));
goodT =  @(X) find(~cellfun(@isempty,X));
bin = 0.01;
filtWidth = 0.02;
pre =  1.0;
post = 1.8;

addpath ../Rajeev_Code
Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/SOMCre/DualModalityCue/';

ExptDates = { '2017-12-18' };

for d = 1:length(ExptDates)
    
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} '/Somcre_mixed.mat'] )
    
    %%
    [Spfc_mix, Smd_mix] = packageData(Z_C1, Smd, Spfc);
    
    %%
    
    ct = 0;
    for i = 1:numel(Smd_mix)
        
        
        clear R1C1 R2C1 spikeRateR1C1 spikeRateR2C1 r1 r2
        % correct =============================================================
        R1C1  = Smd_mix(i).SpikeTimes_R1_nosound;
        R2C1  = Smd_mix(i).SpikeTimes_R2_nosound;
        R1C1  =  R1C1( goodT(R1C1) );
        R2C1  = R2C1( goodT(R2C1) );
        
        if ~isempty(R1C1) && ~isempty(R2C1)
            ct = ct+1;
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            x1 = spikeMatR1C1';
            
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            x2 = spikeMatR2C1';
            
            ind = find( time>=0.25 & time <= 0.35 );
            
            r1 = nanmean( x1(:,ind), 2);
            r2 = nanmean( x2(:,ind), 2);
            
            V(ct) = CD(r1,r2);
            
            X1(ct,:) = spikeRateR1C1;
            X2(ct,:) = spikeRateR2C1;
            
            Sel(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
        end;
    end;
    
    
    %%
    
    ct = 0;
    for i = 1:numel(Smd_mix)
        
        clear R1C1 R2C1 spikeRateR1C1 spikeRateR2C1 r1 r2
        % incorrect ===========================================================
        R1C1  = Smd_mix(i).SpikeTimes_R1_sound;
        R2C1  = Smd_mix(i).SpikeTimes_R2_sound;
        R1C1  =  R1C1( goodT(R1C1) );
        R2C1  = R2C1( goodT(R2C1) );
        
        if ~isempty(R1C1) && ~isempty(R2C1)
            ct = ct+1;
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            x1 = spikeMatR1C1';
            
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            x2 = spikeMatR2C1';
            
            ind = find( time>= 0.25 & time <= 0.35 );
            
            r1 = nanmean( x1(:,ind), 2);
            r2 = nanmean( x2(:,ind), 2);
            
            V_INC(ct) = CD(r1,r2);
            
            X1_INC(ct,:) = spikeRateR1C1;
            X2_INC(ct,:) = spikeRateR2C1;
            
            Sel_INC(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
        end;
        
    end;
    
    indsGood = find( isnan(V)==0 );
    
    Vn = V( indsGood );
    MDCD = Vn./ sum( abs(Vn) );
    
    X1 = X1( indsGood,:);
    X2 = X2( indsGood,:);
    
    MD_proj_x1  = MDCD*X1;
    MD_proj_x2  = MDCD*X2;
    
    
    indsGood = find( isnan(V_INC)==0 );
    
    Vn = V_INC( indsGood );
    MDCD_INC = Vn./ sum( abs(Vn) );
    
    X1_INC = X1_INC( indsGood,:);
    X2_INC = X2_INC( indsGood,:);
    
    MD_proj_x1_INC  = MDCD_INC*X1_INC;
    MD_proj_x2_INC  = MDCD_INC*X2_INC;
    
    
    %%
    figure(1);set(gcf,'color','w');
    subplot(1,2,1)
    plot( time(1:end-1), smooth(MD_proj_x1(1:end-1)),'k'); hold on;
    plot( time(1:end-1), smooth(MD_proj_x2(1:end-1)),'r'); hold on;
    title(['No-Sound'])
    xlim([-1, 1.2]);axis square; box off; set(gca,'tickdir','out','fontsize',16)
    
    subplot(1,2,2)
    plot( time(1:end-1), smooth(MD_proj_x1_INC(1:end-1)),'k'); hold on;
    plot( time(1:end-1), smooth(MD_proj_x2_INC(1:end-1)),'r'); hold on;
    title(['Sound']);
    xlim([-1, 1.2]);axis square; box off; set(gca,'tickdir','out','fontsize',16)
    
    
    figure(2);set(gcf,'color','w');
    plot( time(1:end-1), smooth( MD_proj_x1(1:end-1)- MD_proj_x2(1:end-1)),'g'); hold on;
    plot( time(1:end-1), smooth( MD_proj_x1_INC(1:end-1)-MD_proj_x2_INC(1:end-1)),'b'); hold on;
    xlim([-0.3, 1.2]);axis square; box off; set(gca,'tickdir','out','fontsize',16)
    
    %%
    theta = atand( Sel./Sel_INC );
    R     = sqrt( Sel.^2 + Sel_INC.^2 );
    
    figure(3); set(gcf,'color','w');
    subplot(1,2,1)
    plot( Sel, Sel_INC, 'o', 'markersize', 12, 'markerfacecolor','k'); hold on;
    
    subplot(1,2,2)
    plot( theta, R, 'o', 'markersize', 12, 'markerfacecolor','k'); hold on;

timevec = [-0.25:0.05:1.2];

for t = 1:length(timevec)-1
    clear indsGood V Vn X1 X2
    ct = 0;
    for i = 1:numel(Smd_mix)
        
        
        clear R1C1 R2C1 spikeRateR1C1 spikeRateR2C1 r1 r2
        % sound =============================================================
        R1C1  = Smd_mix(i).SpikeTimes_R1_sound;
        R2C1  = Smd_mix(i).SpikeTimes_R2_sound;
        R1C1  = R1C1( goodT(R1C1) );
        R2C1  = R2C1( goodT(R2C1) );
        
        if ~isempty(R1C1) && ~isempty(R2C1)
            ct = ct+1;
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            x1 = spikeMatR1C1';
            
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            x2 = spikeMatR2C1';
            
            ind = find( time>= timevec(t) & time <= timevec(t+1) );
            
            r1 = nanmean( x1(:,ind), 2);
            r2 = nanmean( x2(:,ind), 2);
            
            V(ct) = CD(r1,r2);
            
            X1(ct,:) = spikeRateR1C1;
            X2(ct,:) = spikeRateR2C1;
            
            Sel(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
        end;
    end;
    
    indsGood = find( isnan(V)==0 );
    
    Vn = V(indsGood);
    MDCD_sound{t} = Vn./ nansum( abs(Vn) );
    
    X1 = X1( indsGood,:);
    X2 = X2( indsGood,:);
    
    MD_sound_proj_x1{t}  = MDCD_sound{t}*X1;
    MD_sound_proj_x2{t}  = MDCD_sound{t}*X2;
    
end;



for t = 1:length(timevec)-1
    clear indsGood V Vn
    ct = 0;
    for i = 1:numel(Spfc_mix)
        
        
        clear R1C1 R2C1 spikeRateR1C1 spikeRateR2C1 r1 r2
        % sound =============================================================
        R1C1  = Spfc_mix(i).SpikeTimes_R1_sound;
        R2C1  = Spfc_mix(i).SpikeTimes_R2_sound;
        R1C1  = R1C1( goodT(R1C1) );
        R2C1  = R2C1( goodT(R2C1) );
        
        if ~isempty(R1C1) && ~isempty(R2C1)
            ct = ct+1;
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            x1 = spikeMatR1C1';
            
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            x2 = spikeMatR2C1';
            
            ind = find( time>= timevec(t) & time <= timevec(t+1) );
            
            r1 = nanmean( x1(:,ind), 2);
            r2 = nanmean( x2(:,ind), 2);
            
            V(ct) = CD(r1,r2);
            
            X1(ct,:) = spikeRateR1C1;
            X2(ct,:) = spikeRateR2C1;
            
            Sel(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
        end;
    end;
    
    indsGood = find( isnan(V)==0 );
    
    Vn = V(indsGood);
    PFCCD_sound{t} = Vn./ nansum( abs(Vn) );
    
    X1 = X1( indsGood,:);
    X2 = X2( indsGood,:);
    
    PFC_sound_proj_x1{t}  = PFCCD_sound{t}*X1;
    PFC_sound_proj_x2{t}  = PFCCD_sound{t}*X2;
    
    
end;


%%

for t = 1:length(timevec)-1
    clear indsGood V Vn
    ct = 0;
    for i = 1:numel(Smd_mix)
        
        
        clear R1C1 R2C1 spikeRateR1C1 spikeRateR2C1 r1 r2
        % sound =============================================================
        R1C1  = Smd_mix(i).SpikeTimes_R1_nosound;
        R2C1  = Smd_mix(i).SpikeTimes_R2_nosound;
        R1C1  = R1C1( goodT(R1C1) );
        R2C1  = R2C1( goodT(R2C1) );
        
        if ~isempty(R1C1) && ~isempty(R2C1)
            ct = ct+1;
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            x1 = spikeMatR1C1';
            
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            x2 = spikeMatR2C1';
            
            ind = find( time>= timevec(t) & time <= timevec(t+1) );
            
            r1 = nanmean( x1(:,ind), 2);
            r2 = nanmean( x2(:,ind), 2);
            
            V(ct) = CD(r1,r2);
            
            X1(ct,:) = spikeRateR1C1;
            X2(ct,:) = spikeRateR2C1;
            
            Sel(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
        end;
    end;
    
    indsGood = find( isnan(V)==0 );
    
    Vn = V(indsGood);
    MDCD_nosound{t} = Vn./ nansum( abs(Vn) );
    
    X1 = X1( indsGood,:);
    X2 = X2( indsGood,:);
    
    MD_nosound_proj_x1{t}  = MDCD_nosound{t}*X1;
    MD_nosound_proj_x2{t}  = MDCD_nosound{t}*X2;
    
end;

for t = 1:length(timevec)-1
    clear indsGood V Vn
    ct = 0;
    for i = 1:numel(Spfc_mix)
        
        
        clear R1C1 R2C1 spikeRateR1C1 spikeRateR2C1 r1 r2
        % sound =============================================================
        R1C1  = Spfc_mix(i).SpikeTimes_R1_nosound;
        R2C1  = Spfc_mix(i).SpikeTimes_R2_nosound;
        R1C1  = R1C1( goodT(R1C1) );
        R2C1  = R2C1( goodT(R2C1) );
        
        if ~isempty(R1C1) && ~isempty(R2C1)
            ct = ct+1;
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            x1 = spikeMatR1C1';
            
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            x2 = spikeMatR2C1';
            
            ind = find( time>= timevec(t) & time <= timevec(t+1) );
            
            r1 = nanmean( x1(:,ind), 2);
            r2 = nanmean( x2(:,ind), 2);
            
            V(ct) = CD(r1,r2);
            
            X1(ct,:) = spikeRateR1C1;
            X2(ct,:) = spikeRateR2C1;
            
            Sel(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
        end;
    end;
    
    indsGood = find( isnan(V)==0 );
    
    Vn = V(indsGood);
    PFCCD_nosound{t} = Vn./ nansum( abs(Vn) );
    
    X1 = X1( indsGood,:);
    X2 = X2( indsGood,:);
    
    PFC_nosound_proj_x1{t}  = PFCCD_nosound{t}*X1;
    PFC_nosound_proj_x2{t}  = PFCCD_nosound{t}*X2;
    
end;

%%

for i = 1:length(timevec)-1
    for j = 1:length(timevec)-1
        foo1 = smooth( MDCD_sound{i} );
        foo2 = smooth( MDCD_sound{j} );
        m = min( length(foo1), length(foo2) );
        MD_CC_Sound(i,j) = corr( foo1(1:m), foo2(1:m) );
    end;
end;

for i = 1:length(timevec)-1
    for j = 1:length(timevec)-1
        foo1 = smooth( MDCD_nosound{i} );
        foo2 = smooth( MDCD_nosound{j} );
        m = min( length(foo1), length(foo2) );
        MD_CC_NoSound(i,j) = corr( foo1(1:m), foo2(1:m) );
    end;
end;

for i = 1:length(timevec)-1
    for j = 1:length(timevec)-1
        foo1 = smooth( PFCCD_sound{i} );
        foo2 = smooth( PFCCD_sound{j} );
        m = min( length(foo1), length(foo2) );
        PFC_CC_Sound(i,j) = corr( foo1(1:m), foo2(1:m) );
    end;
end;

for i = 1:length(timevec)-1
    for j = 1:length(timevec)-1
        foo1 = smooth( PFCCD_nosound{i} );
        foo2 = smooth( PFCCD_nosound{j} );
        m = min( length(foo1), length(foo2) );
        PFC_CC_NoSound(i,j) = corr( foo1(1:m), foo2(1:m) );
    end;
end;

%%
figure(1); set(gcf,'color','w')
subplot(1,2,1);
imagesc(timevec, timevec, MD_CC_Sound); colormap(jet)
axis square; box off; set(gca,'tickdir','out','fontsize',16);
line([-0.2, 1.2],[0,0],'color','k')
line([-0.2, 1.2],[0.25,0.25],'color','k')
line([-0.2, 1.2],[0.35,0.35],'color','k')
line([-0.2, 1.2],[0.85,0.85],'color','k')

line([0,0],[-0.2, 1.2],'color','k')
line([0.25,0.25],[-0.2, 1.2],'color','k')
line([0.35,0.35],[-0.2, 1.2],'color','k')
line([0.85,0.85],[-0.2, 1.2],'color','k')

subplot(1,2,2)
imagesc(timevec, timevec, MD_CC_NoSound); colormap(jet)
axis square; box off; set(gca,'tickdir','out','fontsize',16);
line([-0.2, 1.2],[0,0],'color','k')
line([-0.2, 1.2],[0.25,0.25],'color','k')
line([-0.2, 1.2],[0.35,0.35],'color','k')
line([-0.2, 1.2],[0.85,0.85],'color','k')

line([0,0],[-0.2, 1.2],'color','k')
line([0.25,0.25],[-0.2, 1.2],'color','k')
line([0.35,0.35],[-0.2, 1.2],'color','k')
line([0.85,0.85],[-0.2, 1.2],'color','k')


figure(2); set(gcf,'color','w')
subplot(1,2,1);
imagesc(timevec, timevec, PFC_CC_Sound); colormap(jet)
axis square; box off; set(gca,'tickdir','out','fontsize',16);
line([-0.2, 1.2],[0,0],'color','k')
line([-0.2, 1.2],[0.25,0.25],'color','k')
line([-0.2, 1.2],[0.35,0.35],'color','k')
line([-0.2, 1.2],[0.85,0.85],'color','k')

line([0,0],[-0.2, 1.2],'color','k')
line([0.25,0.25],[-0.2, 1.2],'color','k')
line([0.35,0.35],[-0.2, 1.2],'color','k')
line([0.85,0.85],[-0.2, 1.2],'color','k')

subplot(1,2,2)
imagesc(timevec, timevec, PFC_CC_NoSound); colormap(jet)
axis square; box off; set(gca,'tickdir','out','fontsize',16);
line([-0.2, 1.2],[0,0],'color','k')
line([-0.2, 1.2],[0.25,0.25],'color','k')
line([-0.2, 1.2],[0.35,0.35],'color','k')
line([-0.2, 1.2],[0.85,0.85],'color','k')

line([0,0],[-0.2, 1.2],'color','k')
line([0.25,0.25],[-0.2, 1.2],'color','k')
line([0.35,0.35],[-0.2, 1.2],'color','k')
line([0.85,0.85],[-0.2, 1.2],'color','k')

%%
for i = 1:length(timevec)-1
    PFCNS1_end(i) = max( PFC_nosound_proj_x1{i} );
    PFCNS2_end(i) = min( PFC_nosound_proj_x2{i} );
    
    PFCS1_end(i) = max( PFC_sound_proj_x1{i} );
    PFCS2_end(i) = min( PFC_sound_proj_x2{i} );
    
    MDNS1_end(i) = max( MD_nosound_proj_x1{i} );
    MDNS2_end(i) = min( MD_nosound_proj_x2{i} );
    
    MDS1_end(i) = max( MD_sound_proj_x1{i} );
    MDS2_end(i) = min( MD_sound_proj_x2{i} );
end;
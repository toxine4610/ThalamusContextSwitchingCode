bin = 0.005;
filtWidth = 0.015;
pre =  0.2;
post = 1.0;

CD = @(r1, r2) ( nanmean(r1,1) - nanmean(r2,1) )./sqrt( nanvar(r1,[],1)+nanvar(r2,[],1));


for i = 1:numel(Smd_mix)
    
    R1sound= Smd_mix(i).SpikeTimes_R1_sound;
    R2sound= Smd_mix(i).SpikeTimes_R2_sound;
    
    
    [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1sound, bin, filtWidth, pre, post );
    x1 = spikeMatR1C1';
    
    [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2sound, bin, filtWidth, pre, post );
    x2 = spikeMatR2C1';
    
    ind = find( time>=0 & time <= 0.1 );
    
    r1 = nanmean( x1(1:end,ind), 2);
    r2 = nanmean( x2(1:end-1,ind), 2);
    
    V(i) = CD(r1,r2);
    
    
    X1(i,:) = spikeRateR1C1(1:end-1);
    X2(i,:) = spikeRateR2C1(1:end-1);
end;

indsGood = find( isnan(V)==0 );

Vn = V( indsGood );
MDCD = Vn./ sum( abs(Vn) );

X1 = X1( indsGood,:);
X2 = X2( indsGood,:);

MD_proj_x1  = MDCD*X1;
MD_proj_x2  = MDCD*X2;

figure(1);
plot( time(1:end-1), smooth(MD_proj_x1))
hold on;
plot( time(1:end-1), smooth(MD_proj_x2))

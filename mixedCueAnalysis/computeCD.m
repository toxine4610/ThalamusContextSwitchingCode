
function [MDCD_sound,MD_sound_Sel, MD_sound_proj_x1, MD_sound_proj_x2,...
    MDCD_nosound,MD_nosound_Sel, MD_nosound_proj_x1, MD_nosound_proj_x2,...
    PFCCD_sound,PFC_sound_Sel, PFC_sound_proj_x1, PFC_sound_proj_x2,...
    PFCCD_nosound,PFC_nosound_Sel, PFC_nosound_proj_x1, PFC_nosound_proj_x2, time] = computeCD(Smd_mix, Spfc_mix, timevec)



CD = @(r1, r2) ( nanmean(r1,1) - nanmean(r2,1) )./sqrt( nanvar(r1,[],1)+nanvar(r2,[],1));
goodT =  @(X) find(~cellfun(@isempty,X));

bin = 0.01;
filtWidth = 0.1;
pre =  0.5;
post = 1.8;

% if(~exist('timevec'))
     timevec = linspace(0,0.8,10);
% end


%%
for t = 1:length(timevec)-1
    clear indsGood V Vn X1 X2
    ct = 0;
    for i = 1:numel(Smd_mix)
        
        
        clear R1C1 R2C1 spikeRateR1C1 spikeRateR2C1 r1 r2
         % sound =============================================================
        R1C1  = Smd_mix(i).SpikeTimes_R1_nosound;
        R2C1  = Smd_mix(i).SpikeTimes_R2_nosound;
        R1C1  = R1C1( goodT(R1C1) );
        R2C1  = R2C1( goodT(R2C1) );
        
        R1C2  = Smd_mix(i).SpikeTimes_R1_sound;
        R2C2  = Smd_mix(i).SpikeTimes_R2_sound;
        R1C2  = R1C2( goodT(R1C2) );
        R2C2  = R2C2( goodT(R2C2) );
        
        
        C1 = [R1C1, R2C1];
        C2 = [R1C2, R2C2];
        
        
%         if ~isempty(R1C1) && ~isempty(R2C1)
            ct = ct+1;
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
            x1 = spikeMatR1C1';
            
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
            x2 = spikeMatR2C1';
            
            ind = find( time>= timevec(t) & time <= timevec(t+1) );
            
            r1 = nanmean( x1(:,ind), 2);
            r2 = nanmean( x2(:,ind), 2);
            
            V(ct) = CD(r1,r2);
            
            X1(ct,:) = spikeRateR1C1;
            X2(ct,:) = spikeRateR2C1;
            
            MD_sound_Sel{t}(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
%         end;
    end;
    
    indsGood = find( isnan(V)==0 );
    
    Vn = V(indsGood);
    MDCD_sound{t} = Vn./ nansum( abs(Vn) );
    
    X1 = X1( indsGood,:);
    X2 = X2( indsGood,:);
    
    MD_sound_proj_x1{t}  = MDCD_sound{t}*X1;
    MD_sound_proj_x2{t}  = MDCD_sound{t}*X2;
    
end;



%%
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
        
        R1C2  = Spfc_mix(i).SpikeTimes_R1_sound;
        R2C2  = Spfc_mix(i).SpikeTimes_R2_sound;
        R1C2  = R1C2( goodT(R1C2) );
        R2C2  = R2C2( goodT(R2C2) );
        
        
        C1 = [R1C1, R2C1];
        C2 = [R1C2, R2C2];
        
        
%         if ~isempty(R1C1) && ~isempty(R2C1)
            ct = ct+1;
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
            x1 = spikeMatR1C1';
            
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( C2, bin, filtWidth, pre, post );
            x2 = spikeMatR2C1';
            
            ind = find( time>= timevec(t) & time <= timevec(t+1) );
            
            r1 = nanmean( x1(:,ind), 2);
            r2 = nanmean( x2(:,ind), 2);
            
            V(ct) = CD(r1,r2);
            
            X1(ct,:) = spikeRateR1C1;
            X2(ct,:) = spikeRateR2C1;
            
            PFC_sound_Sel{t}(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
%         end;
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
        R1C1  = Smd_mix(i).SpikeTimes_R1_nosound_inc;
        R2C1  = Smd_mix(i).SpikeTimes_R2_nosound_inc;
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
            
            MD_nosound_Sel{t}(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
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
        R1C1  = Spfc_mix(i).SpikeTimes_R1_nosound_inc;
        R2C1  = Spfc_mix(i).SpikeTimes_R2_nosound_inc;
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
            
            PFC_nosound_Sel{t}(ct) = ( nanmean(r1)-nanmean(r2) )./( nanmean(r1)+nanmean(r2) );
            
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

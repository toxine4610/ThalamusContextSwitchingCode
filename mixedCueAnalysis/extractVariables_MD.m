


function Out = extractVariables_MD( SMD_mix )

%%
CMI   = @(X,Y) (X-Y)./(X+Y);
goodT =  @(X) find(~cellfun(@isempty,X));
bin = 0.01;
filtWidth = 0.05;
pre =  0.2;
post = 1.8;

%%
ct = 0;
for i = 1:numel(SMD_mix)
%     try
        clear C1 C2 R1C1 R2C1 R1C2 R2C2 spikeRateR1C1 spikeRateR2C1 spikeRateR1C2 spikeRateR2C2
        % correct =============================================================
       
        if ~isempty(SMD_mix(i).SpikeTimes_R1_sound)
            R1C1  = SMD_mix(i).SpikeTimes_R1_sound;
            R2C1  = SMD_mix(i).SpikeTimes_R2_sound;
        elseif ~isempty(SMD_mix(i).SpikeTimes_R1_nosound)
            R1C1  = SMD_mix(i).SpikeTimes_R1_nosound;
            R2C1  = SMD_mix(i).SpikeTimes_R2_nosound;
        end
        
        R1C1  =  R1C1( goodT(R1C1) );
        R2C1  =  R2C1( goodT(R2C1) );
        
        C1 = [R1C1, R2C1];
        
        [spikeRateC1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( C1, bin, filtWidth, pre, post );
        
        MDFRTS_C1(i,:) = spikeRateC1;
        
  
        indPCD  = find( time >= 0 & time <= 0.1 );
        indBase = find( time >= -0.25 & time  <=0 );
        indDelay = find( time >=0.1 & time <=0.50 );
        
        MD_muFR_PCD_C1(i) = nanmean(spikeRateC1(indPCD));
        
        MD_muFR_CFD_C1(i) = nanmean(spikeRateC1(indDelay));
        
        MD_muFR_BAS_C1(i) = nanmean(spikeRateC1(indBase));
        
        MD_CMI_C1_we(i)  =  CMI( MD_muFR_PCD_C1(i), MD_muFR_CFD_C1(i));
        %
        if size(R1C1,2) ~= 0
            [spikeRateR1C1,errorBounds,spikeMatR1C1,time] = estimateSpikeRates( R1C1, bin, filtWidth, pre, post );
            SM11 =  spikeMatR1C1';
            
        else
            spikeRateR1C1 = NaN(1,201);
            SM11 = NaN(22, 201) ;
            
        end;
        
        if length(spikeRateR1C1) == 1
            spikeRateR1C1 = NaN(1,length(time));
        end;

        MD_muFR_PCD_R1C1(i) = nanmean(spikeRateR1C1(indPCD));
        
        MD_muFR_CFD_R1C1(i) = nanmean(spikeRateR1C1(indDelay));
        
        MD_muFR_BAS_R1C1(i) = nanmean(spikeRateR1C1(indBase));
        
        
        if size(R2C1,2) ~= 0
            [spikeRateR2C1,errorBounds,spikeMatR2C1,time] = estimateSpikeRates( R2C1, bin, filtWidth, pre, post );
            SM21 =  spikeMatR2C1';
        else
            spikeRateR2C1 = NaN(1,201);
            SM21 = NaN(22, 201);
        end;
        
        
        if length(spikeRateR2C1) == 1
            spikeRateR2C1 = NaN(1,length(time));
        end;
        
        
        MD_muFR_PCD_R2C1(i) = nanmean(spikeRateR2C1(indPCD));
        
        MD_muFR_CFD_R2C1(i) = nanmean(spikeRateR2C1(indDelay));
        
        MD_muFR_BAS_R2C1(i) = nanmean(spikeRateR2C1(indBase));
        
        MD_Sel1(i) =  ( MD_muFR_CFD_R1C1(i) -  MD_muFR_PCD_R1C1(i));
        MD_Sel2(i) =  (MD_muFR_CFD_R2C1(i) -  MD_muFR_PCD_R2C1(i));
        
        RMD(i)   = sqrt( MD_Sel1(i).^2  +  MD_Sel2(i).^2 );
        ThetaMD(i) = atan( MD_Sel1(i)./ MD_Sel2(i) );
        
        
        MD_dp_C1(i,:) = (nanmean( SM11 ) - nanmean( SM21 ))./ sqrt( 0.5.*( var( SM21,[],1) + var( SM11,[], 1) ) );
       
        
        
%     catch err
%     end
end;

%%
Out.MD_dp_C1 = MD_dp_C1;
Out.RMD      = RMD;
Out.ThetaMD  = ThetaMD;
Out.MD_Sel1  = MD_Sel1;
Out.MD_Sel2 =  MD_Sel2;
Out.MD_muFR_PCD_R2C1 =  MD_muFR_PCD_R2C1;
Out.MD_muFR_CFD_R2C1 =  MD_muFR_CFD_R2C1;
Out.MD_muFR_BAS_R2C1 =  MD_muFR_BAS_R2C1;
Out.MD_muFR_PCD_R1C1 =  MD_muFR_PCD_R1C1;
Out.MD_muFR_CFD_R1C1 =  MD_muFR_CFD_R1C1;
Out.MD_muFR_BAS_R1C1 =  MD_muFR_BAS_R1C1;
Out.MD_muFR_PCD_C1   = MD_muFR_PCD_C1;
Out.MD_muFR_CFD_C1   = MD_muFR_CFD_C1;
Out.MD_muFR_BAS_C1   = MD_muFR_BAS_C1;

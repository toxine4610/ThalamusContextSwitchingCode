
function Z = getZ( filename, MouseID );

[ initT, trial, laser, LEDTBAud, right_trials, left_trials, stim_lengths, catchh, audioStim, visStim, audioStimDistract, visStimDistract, corr, inc, ts, cueL, cueR, ITI, LEDstate] = parseText( filename, MouseID);


offSet = 0 ;
negOff = 0 ;

 
    if offset > 0
        if negOff <= 0
            visStims = intersect(visStimDistract(find(visStimDistract<=(length(initT)+negOff))),visStimDistract(find(visStimDistract >(-offset+1))))
            audStims = intersect(audStimDistract(find(audStimDistract<=(length(initT)+negOff))),audStimDistract(find(audStimDistract > (-offset+1))))
            noStims = intersect(visStim(find(visStim<=(length(initT)+negOff))),visStim(find(visStim >(-offset+1))))
            
            ind1 = find(isrepeat_same == 1);
            repssame = intersect(ind1(find(ind1<=(length(initT)+negOff))), ind1(find(ind1 >(-offset+1))));
            ind2 = find(isrepeat_different == 1);
            repdiff = intersect(ind2(find(ind2<=(length(initT)+negOff))), ind2(find(ind2 >(-offset+1))));
            
            ind1 = find(isdim == 1);
            repdim = intersect(ind1(find(ind1<=(length(initT)+negOff))), ind1(find(ind1 >(-offset+1))))
            
            if length(las) > 0
                las = intersect(las(find(las<=(length(initT)+negOff))),las(find(las > (-offset+1))));
            end
%             Linit.start_time = L2.start_time((offset+1):length(L2.start_time)) - initT((1):(length(initT)+negOff))'/1000
%             testDiff = L2.start_time((offset+1):length(L2.start_time)) - initT(1:(length(initT)+negOff))'/1000
        else
            visStims = visStimDistract;
            audStims = audStimDistract;
            noStims = visStim;
            ind1 = find(isrepeat_same == 1);
            repssame = ind1;
            ind2 = find(isrepeat_different == 1);
            repdiff = ind2;
            ind1 = find(isdim == 1);
            repdim = ind1;
            
            Linit.start_time = L2.start_time((offset+1):(length(L2.start_time)-negOff)) - initT(:)'/1000
            %         Linit.start_time = L2.start_time( 119 : end );
            testDiff = L2.start_time((offset+1):(length(L2.start_time)-negOff)) - initT(:)'/1000
        end
        fVisCorr = (intersect(visStims,corrV));
        fAudCorr = (intersect(audStims,corrV));
        fVisInc = (intersect(visStims,inc));
        fAudInc = (intersect(audStims,inc));
        fRight = [intersect(right_trials,visStims) intersect(right_trials,audStims)]
        
        fLongLatAud = (intersect(audStims,(find(latency > prctile(latency,60)))));
        fLongLatVis = (intersect(visStims,(find(latency > prctile(latency,60)))));
        
        fShortLatAud = (intersect(audStims,(find(latency < prctile(latency,60)))));
        fShortLatVis = (intersect(visStims,(find(latency < prctile(latency,60)))));
        
        allStims = [visStims audStims noStims];
        latency = latency(sortrows(allStims'));
        ts = ts(sortrows(allStims'));
        
    elseif offset < 0
        if negOff < 0
            
            
            visStims = intersect(visStimDistract(find(visStimDistract<=(length(initT)+negOff))),visStimDistract(find(visStimDistract >(-offset))))
            audStims = intersect(audStimDistract(find(audStimDistract<=(length(initT)+negOff))),audStimDistract(find(audStimDistract > (-offset))))
            noStims = intersect(visStim(find(visStim<=(length(initT)+negOff+1))),visStim(find(visStim > (-offset+1)))) %catch trial
            Linit.start_time = L2.start_time(1:length(L2.start_time))' - initT((-offset+1):(length(initT)+negOff))'/1000
            testDiff = L2.start_time(1:length(L2.start_time))' - initT((-offset+1):(length(initT)+negOff))'/1000
            
            ind1 = find(isrepeat_same == 1);
            repssame = intersect(ind1(find(ind1<=(length(initT)+negOff))),ind1(find(ind1 >(-offset))))
            ind2 = find(isrepeat_different == 1);
            repdiff = intersect(ind2(find(ind2<=(length(initT)+negOff))),ind2(find(ind2 >(-offset))))
            
            ind1 = find(isdim == 1);
            repdim = intersect(ind1(find(ind1<=(length(initT)+negOff))),ind1(find(ind1 >(-offset))))
            
        else
            
            visStims = visStimDistract(find(visStimDistract >=(-offset+1)));
            audStims = audStimDistract(find(audStimDistract >=(-offset+1)));
            noStims = visStim(find(visStim >(-offset+1))); %catch trial
            Linit.start_time = L2.start_time(1:(length(L2.start_time)-negOff-1)) - initT((-offset+1):(length(initT)-1))'/1000;
            testDiff = L2.start_time(1:(length(L2.start_time)-negOff-1)) - initT((-offset+1):(length(initT)-1))'/1000;
            
            ind1 = find(isrepeat_same == 1);
            repssame = ind1( find(ind1 >=(-offset+1)) );
            ind2 = find(isrepeat_different == 1);
            repdiff = ind2( find(ind2 >=(-offset+1)) );
            
            ind1 = find(isdim == 1);
            repdim = ind1( find(ind1 >=(-offset+1)) );
        end
        
        
        allStims = [visStims audStims noStims];
        allStims = sortrows(allStims')';
        allStims = allStims(1:length(Linit.start_time));
        
        latency = latency(sortrows(allStims'));
        
        corrV = intersect(corrV,allStims);
        inc = intersect(inc,allStims);
        if length(las) > 0
            las = intersect(las,allStims);
        end
        ts = ts(sortrows((allStims)'));
        allStims = (allStims);
        visStims = (visStims);
        audStims = (audStims);
        
        
        fVisCorr = (intersect(visStims,corrV));
        fAudCorr = (intersect(audStims,corrV));
        fVisInc = (intersect(visStims,inc));
        fAudInc = (intersect(audStims,inc));
        fRight = [intersect(right_trials,visStims) intersect(right_trials,audStims)];
        
        allStims = [fVisCorr,fAudCorr,fVisInc,fAudInc,noStims];
        allStims = sortrows(allStims')';
        
        fLongLatAud = (intersect(audStims,(find(latency > prctile(latency,75)))));
        fLongLatVis = (intersect(visStims,(find(latency > prctile(latency,75)))));
        
        fShortLatAud = (intersect(audStims,(find(latency < prctile(latency,75)))));
        fShortLatVis = (intersect(visStims,(find(latency < prctile(latency,75)))));
    else % offeset == 0
        if negOff < 0
            
                        initT = initT(1:end-1);

            visStims = intersect(visStimDistract(find(visStimDistract<=(length(initT)+negOff+1))),visStimDistract(find(visStimDistract >(offset))))
            audStims = intersect(audStimDistract(find(audStimDistract<=(length(initT)+negOff+1))),audStimDistract(find(audStimDistract > (offset))))
            noStims = intersect(visStim(find(visStim<=(length(initT)+negOff+1))),visStim(find(visStim > (offset)))) %catch trial
            Linit.start_time = L2.start_time((offset+1):length(L2.start_time)) - initT((1):(length(initT)+negOff))'/1000
            testDiff = L2.start_time((offset+1):length(L2.start_time)) - initT(1:(length(initT)+negOff))'/1000
            ind1 = find(isrepeat_same == 1);
            repssame = intersect(ind1(find(ind1<=(length(initT)+negOff))),ind1(find(ind1 >(-offset+1))));
            ind2 = find(isrepeat_different == 1);
            repdiff = intersect(ind2(find(ind2<=(length(initT)+negOff))),ind2(find(ind2 >(-offset+1))));
            
            ind1 = find(isdim == 1);
            repdim = intersect(ind1(find(ind1<=(length(initT)+negOff))),ind1(find(ind1 >(-offset+1))));
            
        else
            visStims = visStimDistract(find(visStimDistract >=(offset)));
            audStims = audStimDistract(find(audStimDistract >=(offset)));
            noStims = visStim(find(visStim >(offset))); %catch trial
            Linit.start_time = L2.start_time(1:(length(L2.start_time)-negOff));% - initT((-offset+1):(length(initT)))'/1000;
            %          Linit.start_time = L2.start_time( 119: end );
            %          testDiff =   L2.start_time( 119 : end );
            testDiff = L2.start_time(1:(length(L2.start_time)-negOff));% - initT((-offset+1):(length(initT)))'/1000;
            
            ind1 = find(isrepeat_same == 1);
            repssame = ind1(find( ind1 >=(offset)));
            ind2 = find(isrepeat_different == 1);
            repdiff = ind2(find( ind2 >=(offset)));
            
            ind1 = find(isdim == 1);
            repdim = ind1(find( ind1 >=(offset)));
            
        end
        
        
        allStims = [visStims audStims noStims];
        allStims = sortrows(allStims')';
        allStims = allStims(1:length(Linit.start_time));
        
        %     isRep_same = isrepeat_Same(1:length(Linit.start_time));
        %     isRep_diff = isrepeat_Conflict(1:length(Linit.start_time));
        
        latency = latency(sortrows(allStims'));
        
        corrV = intersect(corrV,allStims) + offset;
        inc = intersect(inc,allStims)+offset;
        if length(las) > 0
            las = intersect(las,allStims)+offset;
        end
        ts = ts(sortrows((allStims)'));
        allStims = (allStims+offset);
        visStims = (visStims+offset);
        audStims = (audStims+offset);
        
        
        fVisCorr = (intersect(visStims,corrV));
        fAudCorr = (intersect(audStims,corrV));
        fVisInc = (intersect(visStims,inc));
        fAudInc = (intersect(audStims,inc));
        fRight = [intersect(right_trials,visStims) intersect(right_trials,audStims)];
        
        allStims = [fVisCorr,fAudCorr,fVisInc,fAudInc,noStims];
        allStims = sortrows(allStims')';
        
        fLongLatAud = (intersect(audStims,(find(latency > prctile(latency,75)))))+offset;
        fLongLatVis = (intersect(visStims,(find(latency > prctile(latency,75)))))+offset;
        
        fShortLatAud = (intersect(audStims,(find(latency < prctile(latency,75)))))+offset;
        fShortLatVis = (intersect(visStims,(find(latency < prctile(latency,75)))))+offset;
    end
    
    fVisCorr = fVisCorr(find(fVisCorr < length(allStims)));% original error from here offset?
    fAudCorr = fAudCorr(find(fAudCorr < length(allStims)));
    
    fRight = fRight(find(fRight < length(allStims)));
    
    outCome = zeros(length(allStims),1);
    outCome(fVisCorr) = 1;%+offset
    outCome(fAudCorr) = 1;
    %outCome(noStims-1) = NaN;
    trialType = zeros(length(allStims),1);
    trialType(visStims) = 1;
    trialType(audStims) = 2;
    trialType = trialType(1:length(allStims));
    latency = latency(1:length(allStims));
    latencysizeTest = size(latency)
    if latencysizeTest(1) < latencysizeTest(2)
        latency = latency'
    end
    % Z = [outCome trialType latency'];
    Z = [outCome trialType latency];
    
    Z(:,5) = 0;
    Z(fRight,5) = 1;
    %Z(noStims-1,5) = NaN;
    
    if length(las) > 0
        lasL = outCome*0;
        lasL(las) = 1;
        lasL = lasL(1:length(allStims));
        Z = [Z lasL]
    else
        Z(:,6) = 0
    end
    
    tfdir=isdir(figpath);
    if tfdir==0
        mkdir(figpath);
    end
    cd(figpath);
    figure;subplot(2,1,1);plot(diff(testDiff))
    subplot(2,1,2);plot(diff(ts))
    saveas(gcf,'testDiff-tsdiff.fig','fig');
    % close all;
    
    hold_time = initT/1000;
    
    length(fVisInc)
    length(fVisCorr)
    
    LnewHold = Linit;
    
    Z(:,4) = 0;
    
    trial = trial(allStims);
    
    Z(:,7) = trial;
    w=gausswin(7);
    w=w/sum(w);
    Z(:,8)=conv(Z(:,1),w,'same');
    
    repType = zeros(length(allStims),1);
    repType(repssame) = 1;
    repType(repdiff)  = 2;
    repType(repdim)   = 3;
    repType(setdiff( (1:length(allStims)), [repdiff, repssame, repdim] )) = 0;
    repType = repType(1:length(allStims));
    
    Z(:,9) = repType;
    
    % Z(:,9) = isrepeat_same(1:length(Linit.start_time));
    % Z(:,10) = isrepeat_different(1:length(Linit.start_time));
    
    
    Zholder = Z;
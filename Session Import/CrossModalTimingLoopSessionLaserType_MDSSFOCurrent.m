clear;
path= 'Y:\CheetahData\WT24\2017-11-03_11-37-42';
files=dir(path);
cd(path);
extt = 1;

% z=0;
% for j=3:length(files)
%     filename=files(j).name;
%     if filename(end-5:end)=='er.mat'
%         z=z+1
%         TTn=filename(1:find(filename=='c')-1);
%         TT=load(filename);
%         unt(z)=length(TT.(TTn).labels);
%     end
% end
%%
clear;
for sn = [36]%0
    for  cn = [1]%
        for tn=[1:16]
            
            extt = 0
        if cn==1
            laserinit=1;%amp4 31
        else
            laserinit=1;%amp4 21
        end
%%
% clear all
% sn=29;
% cn=2;
% laserinit=1;

% initfilePFCWT19MDSSFO2;
% initfilePFCMD128_Olivier;
% initfilePFCamp4_Olivier
% intifilePFCMD_Guang_WT24;

intifilePFCMD_Guang;
session_num =  sn;
%all_tt_nums = [2 9 10 11 13 14 15 18 19 23 26 28 31 32 36 37 38 43 44];
Fpath=Se.folder{session_num};
cd(Se.folder{session_num})

%%
% names = dir('TT*')
% newNames = names;
% 
% for n = 1:length(names)
%     newNames(n).name = ['Sc' newNames(n).name(3:end)]
%     dos(['rename "' names(n).name '" "' newNames(n).name '"']);
% end
% 
% 
% names = dir('ST*')
% newNames = names;
% 
% for n = 1:length(names)
%     newNames(n).name = ['Sc' newNames(n).name(3:end)]
%     dos(['rename "' names(n).name '" "' newNames(n).name '"']);
% end

% names = dir('*.Nst')
% newNames = names;
% 
% for n = 1:length(names)
%     newNames(n).name(find(names(n).name == 'T')) = 'c';
% %     newNames(n).name(find(names(n).name == 'n')) = 'N';
% %     newNames(n).name(find(names(n).name == 'N')) = 'n';
% %     newNames(n).name(find(names(n).name == 's')) = 't';
%     dos(['rename "' names(n).name '" "' newNames(n).name '"']);
% end
% 
% names = dir('*.ncs')
% newNames = names;
% 
% for n = 1:length(names)
%     newNames(n).name(find(names(n).name == 'n')) = 'N';
%     dos(['rename "' names(n).name '" "' newNames(n).name '"']);
% end

% names = dir('*.cut')
% newNames = names;
% 
% for n = 1:length(names)
%     newNames(n).name(find(names(n).name == 'T')) = 'c';
%     newNames(n).name(find(names(n).name == 'n')) = 'N';
%     newNames(n).name(1) = 'S';
%     dos(['rename "' names(n).name '" "' newNames(n).name '"']);
% end
% names = dir('*.ntt')
% newNames = names;
% 
% for n = 1:length(names)
%     newNames(n).name(1) = 'S';
%     newNames(n).name(2) = 'c';
%     newNames(n).name(find(names(n).name == 'N')) = 'n';
% %    newNames(n).name(find(names(n).name == 't')) = 'st';
%     dos(['rename "' names(n).name '" "' newNames(n).name '"']);
% end


%%


% names = dir('*.cut');
% 
% all_tt_nums = [];
% 
% for n = 1:length(names)
%     nameSU = names(n).name;
%     nameSU(strfind(nameSU,'.cut'):length(nameSU)) = [];
%     nameSU(1:2) = [];
%     all_tt_nums = [all_tt_nums str2num(nameSU)];
% end
% 
% all_tt_nums = sort(all_tt_nums);
% 
% eventData = read_cheetah_data([Se.folder{session_num} '\Events.Nev'])
% % 
%  i = 0;
% % 
% % TTLstim1 = 16;
% % 
% for n = 1:(length(all_tt_nums))
%     curr_tt_num = all_tt_nums(n);
%     Sc = spikes(Se,session_num,curr_tt_num);
%     
%     obj.folder = [Fpath '\Sc' num2str(curr_tt_num) '.cut'];%TT/Sc
%     cluster_info = load([obj.folder '']);
%     maxun=max(cluster_info);
%     
%     m = 0;
%     for mn=[1:maxun]
%         cl_holder = cluster(Sc,mn);
%         is_full = max(size(cl_holder.timestamp));
%         if is_full > 1
%             i = i+1;
%             m = m+1;
%             eval(sprintf('cl%d = cluster(Sc,mn)', i));close all
%             num_seq(i,1:2) = [curr_tt_num mn];
%         end
%     end
%     Sc_unit_count(n,1) = curr_tt_num;
%     Sc_unit_count(n,2) = m;
%     
%     
%     is_cluster = 1;
%     m = 0;
%     while is_cluster == 1
%         m = m+1;
%         cl_holder = cluster(Sc,m);
%         is_full = max(size(cl_holder.timestamp));
%         if is_full > 1
%             i = i+1;
%             eval(sprintf('cl%d = cluster(Sc,m)', i));close all
%             num_seq(i,1:2) = [curr_tt_num m];
%         else
%             m = m-1;
%             is_cluster = 0;
%         end
%     end
%     Sc_unit_count(n,1) = curr_tt_num;
%     Sc_unit_count(n,2) = m;
% end
% clear cl_holder i curr_tt_num is_cluster is_full m n mn

%%
% all_tt_nums = [];
% 
% 
% cd(Se.folder{session_num})
% mkdir('AnalyzedFiles')
% names = dir('*.cut');
% 
% for n = 1:length(names)
%     nameSU = names(n).name;
%     nameSU(strfind(nameSU,'.cut'):length(nameSU)) = [];
%     nameSU(1:2) = [];
%     all_tt_nums = [all_tt_nums str2num(nameSU)];
% end
% 
% all_tt_nums = sort(all_tt_nums); %Automatic selection of electrode sets to plot



%%
unitCount = 1;

eventData = read_cheetah_data([Se.folder{session_num} '\Events.nev'])

autoClusterSpikeFiles = dir('TT*cluster.mat');

for n = 1:length(autoClusterSpikeFiles)
    TTname = autoClusterSpikeFiles(n).name;
    numTT(n) = str2num(TTname(3:(length(TTname)-11)));
end

[TTNums,TTOrder] = sort(numTT);

for i = tn % 1:length(TTNums)
    fname = ['TT' num2str(TTNums(i)) '.ntt'];
    tetrodeSpikeTimes = read_cheetah_data(fname);
    load([autoClusterSpikeFiles(TTOrder(i)).name])
    eval(['tetrodeSet = TT' num2str(TTNums(i)) '.labels;'])
    Sc = spikesAuto(Se,session_num,TTNums(i));
    for ii = 1:length(tetrodeSet)
        if max(tetrodeSet{ii}) < length(tetrodeSpikeTimes.ts)
            cl_holder = clusterAuto(Sc,ii);
            unitSpikeTimes = tetrodeSpikeTimes.ts(tetrodeSet{ii});
            eval(['unitAC' num2str(unitCount) ' = unitSpikeTimes;'])
            eval(['cl' num2str(unitCount) ' = cl_holder;'])
            autoClusterSpikeTimes{i}{ii} = unitSpikeTimes;
            num_seq(unitCount,:) = [(TTNums(i)) ii];
        else
            eval(['tetrodeSpikeTimesAuto = TT' num2str(TTNums(i)) '.timestamps;'])
            cl_holder = clusterAuto(Sc,ii);
            unitSpikeTimes = tetrodeSpikeTimesAuto(tetrodeSet{ii});
            eval(['unitAC' num2str(unitCount) ' = unitSpikeTimes;'])
            eval(['cl' num2str(unitCount) ' = cl_holder;'])
            autoClusterSpikeTimes{i}{ii} = unitSpikeTimes;
            num_seq(unitCount,:) = [(TTNums(i)) ii];
        end          
            unitCount = unitCount+1;
    end
end


%%
[TTLc1, TTLc2, TTLn] = TTLlist(sn);
if cn==1% C1 is the first recorded context
 TTLstim1 = TTLc1(1);%laser (1 or 128)
TTLbehave =  2; %stim(vis/aud) (4 or 64) (4 context1 4-poke; 64 context2 3-poke)
else
 TTLstim1 = TTLc2(1);%laser (1 or 128)
TTLbehave = 2;%stim(vis/aud) (4 or 64) (4 context1 4-poke; 64 context2 3-poke) %reverse pin num if wrong
end    
TTLbehave2 = 193;
% behaveEndTrial = 383;

eventData = read_cheetah_data([Se.folder{session_num} '\Events.nev']);
% eventData.TTLval(find(eventData.TTLval(1:452)==64))=0;
eventDataholder = eventData;
% eventData.ts = eventData.ts(1:behaveEndTrial);
% eventData.TTLval = eventData.TTLval(1:behaveEndTrial);
Lstim.start_time = eventData.ts(find(eventData.TTLval == TTLstim1))';
Lstim.pulse_width = diff(eventData.ts);
Lstim.pulse_width = Lstim.pulse_width(find(eventData.TTLval == TTLstim1))';
L.start_time = eventData.ts(find((eventData.TTLval == TTLbehave)+(eventData.TTLval == TTLbehave2)))';
% L.start_time = L.start_time(8:end);
%L.start_time = eventData.ts(find(eventData.TTLval == TTLbehave))';

% fileStruct = dir('*.txt');
% filename = fileStruct.name;
behaviorfilename = behaviorfiles(sn);

% if cn==1
%     [initT, trial, las, LEDTBAud, right_trials, left_trials, stim_lengths, catchh, audioStim, visStim, audStimDistract, visStimDistract, corrV, inc, ts, cueL, cueR, ITI, LEDstate] = parseEngineCrossModalPFCMod(behaviorfilename.C1);
%     initT = initT(1:length(ts))';
% else
%     [initT, trial, las, LEDTBAud, right_trials, left_trials, stim_lengths, catchh, audioStim, visStim, audStimDistract, visStimDistract, corrV, inc, ts, cueL, cueR, ITI, LEDstate] = parseEngineCrossModalPFCMod(behaviorfilename.C2);
%     initT = initT(1:length(ts))';
% end

[initT, trial, las, LEDTBAud, right_trials, left_trials, stim_lengths, catchh, audioStim, visStim, audStimDistract, visStimDistract, corrV, inc, ts, cueL, cueR, ITI1, LEDstate] = parseEngineCrossModalPFCMod(behaviorfilename.C1);
initT = initT(1:length(ts))';
[initT, trial, las, LEDTBAud, right_trials, left_trials, stim_lengths, catchh, audioStim, visStim, audStimDistract, visStimDistract, corrV, inc, ts, cueL, cueR, ITI2, LEDstate] = parseEngineCrossModalPFCMod(behaviorfilename.C2);
initT = initT(1:length(ts))';

if abs(length(L.start_time)-length(ITI1))<=abs(length(L.start_time)-length(ITI2))
    [initT, trial, las, LEDTBAud, right_trials, left_trials, stim_lengths, catchh, audioStim, visStim, audStimDistract, visStimDistract, corrV, inc, ts, cueL, cueR, ITI, LEDstate] = parseEngineCrossModalPFCMod(behaviorfilename.C1);
    initT = initT(1:length(ts))';
else
    [initT, trial, las, LEDTBAud, right_trials, left_trials, stim_lengths, catchh, audioStim, visStim, audStimDistract, visStimDistract, corrV, inc, ts, cueL, cueR, ITI, LEDstate] = parseEngineCrossModalPFCMod(behaviorfilename.C2);
    initT = initT(1:length(ts))';
end
abs(length(L.start_time)-length(ITI1))
abs(length(L.start_time)-length(ITI2))
% trial=[trial trial(end)];
% clear ITI1 ITI2

% [isrepeat_same, isrepeat_different, norepeat] = parseLED(LEDstate);
[isrepeat_same, isrepeat_different, norepeat] = findRepeats(LEDstate);

%%  Visual vs Auditory Behavior Aligner
offsetCorrect = 0
if offsetCorrect
    initT = initT/2;
end
latency = initT;
% 
L.start_time = L.start_time(find(diff(L.start_time) >8));
latency = initT;

sInitTest = size(initT);
if sInitTest(1) < sInitTest(2)
    initT = initT';
end
%L = Lholder;

[a,b,offset] = alignsignals(round(diff(ts/1000)),round(diff(L.start_time)));
offset

% if offset == 0
%     if length(L.start_time) > length(ts)
%          compDiff = [(L.start_time(1)-ts(1)/1000) (L.start_time(2)-ts(2)/1000)]
%         for ii = 1:length(ts)
%             for iii = 1:length(L.start_time)
%                 diffVal(iii) = (mean(compDiff)-(L.start_time(iii)-ts(ii)/1000));
%             end
%             LholdnewTest(ii) =L.start_time(min(find(abs(diffVal) == min(abs(diffVal)))));
%             L.start_time(ii) = L.start_time(min(find(abs(diffVal) == min(abs(diffVal)))));
%         end
%     end
% end
%L.start_time = L.start_time(find(diff(L.start_time) >10))    
[a,b,negOff] = alignsignals(round(diff(wrev(ts)/1000)),round(diff(wrev(L.start_time))));
negOff
offset


if offset > 0
    if negOff <= 0
        visStims = intersect(visStimDistract(find(visStimDistract<=(length(initT)+negOff))),visStimDistract(find(visStimDistract >(-offset+1))))
        audStims = intersect(audStimDistract(find(audStimDistract<=(length(initT)+negOff))),audStimDistract(find(audStimDistract > (-offset+1))))
        noStims = intersect(visStim(find(visStim<=(length(initT)+negOff))),visStim(find(visStim >(-offset+1))))
        if length(las) > 0
            las = intersect(las(find(las<=(length(initT)+negOff))),las(find(las > (-offset+1))));
        end
        Linit.start_time = L.start_time((offset+1):length(L.start_time)) - initT((1):(length(initT)+negOff))'/1000
        testDiff = L.start_time((offset+1):length(L.start_time)) - initT(1:(length(initT)+negOff))'/1000
    else
        visStims = visStimDistract;
        audStims = audStimDistract;
        noStims = visStim;
        Linit.start_time = L.start_time((offset+1):(length(L.start_time)-negOff)) - initT(:)'/1000
%         Linit.start_time = L.start_time( 119 : end );
        testDiff = L.start_time((offset+1):(length(L.start_time)-negOff)) - initT(:)'/1000
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
        Linit.start_time = L.start_time(1:length(L.start_time))' - initT((-offset+1):(length(initT)+negOff))'/1000
        testDiff = L.start_time(1:length(L.start_time))' - initT((-offset+1):(length(initT)+negOff))'/1000
    else
        visStims = visStimDistract(find(visStimDistract >=(-offset+1)));
        audStims = audStimDistract(find(audStimDistract >=(-offset+1)));
        noStims = visStim(find(visStim >(-offset+1))); %catch trial
        Linit.start_time = L.start_time(1:(length(L.start_time)-negOff-1)) - initT((-offset+1):(length(initT)-1))'/1000;
        testDiff = L.start_time(1:(length(L.start_time)-negOff-1)) - initT((-offset+1):(length(initT)-1))'/1000;
        
        ind1 = find(isrepeat_same == 1);
        repssame = ind1( find(ind1 >=(-offset+1)) );
        ind2 = find(isrepeat_different == 1);
        repdiff = ind2( find(ind2 >=(-offset+1)) );
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
        visStims = intersect(visStimDistract(find(visStimDistract<=(length(initT)+negOff+1))),visStimDistract(find(visStimDistract >(offset))))
        audStims = intersect(audStimDistract(find(audStimDistract<=(length(initT)+negOff+1))),audStimDistract(find(audStimDistract > (offset))))
        noStims = intersect(visStim(find(visStim<=(length(initT)+negOff+1))),visStim(find(visStim > (offset)))) %catch trial
        Linit.start_time = L.start_time((offset+1):length(L.start_time)) - initT((1):(length(initT)+negOff))'/1000
        testDiff = L.start_time((offset+1):length(L.start_time)) - initT(1:(length(initT)+negOff))'/1000
        ind1 = find(isrepeat_same == 1);
        repssame = intersect(ind1(find(ind1<=(length(initT)+negOff))),ind1(find(ind1 >(-offset+1))));
        ind2 = find(isrepeat_different == 1);
        repdiff = intersect(ind2(find(ind2<=(length(initT)+negOff))),ind2(find(ind2 >(-offset+1))));

    else
        visStims = visStimDistract(find(visStimDistract >=(offset)));
        audStims = audStimDistract(find(audStimDistract >=(offset)));
        noStims = visStim(find(visStim >(offset))); %catch trial
        Linit.start_time = L.start_time(1:(length(L.start_time)-negOff)) - initT((-offset+1):(length(initT)))'/1000;
%          Linit.start_time = L.start_time( 119: end );
%          testDiff =   L.start_time( 119 : end );
        testDiff = L.start_time(1:(length(L.start_time)-negOff)) - initT((-offset+1):(length(initT)))'/1000;
        
        ind1 = find(isrepeat_same == 1);
        repssame = ind1(find( ind1 >=(offset)));
        ind2 = find(isrepeat_different == 1);
        repdiff = ind2(find( ind2 >=(offset)));
        
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

if cn==1
    figpath=strcat(Fpath,'\GlobalpeakC1');
else
    figpath=strcat(Fpath,'\GlobalpeakC2');
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

Z(:,4) = Linit.start_time(1:length(allStims));

trial = trial(allStims);

Z(:,7) = trial;
w=gausswin(7);
w=w/sum(w);
Z(:,8)=conv(Z(:,1),w,'same');

repType = zeros(length(allStims),1);
repType(repssame) = 1;
repType(repdiff) = 2;
repType(setdiff( (1:length(allStims)), [repdiff, repssame] )) = 0;
repType = repType(1:length(allStims));

Z(:,9) = repType;

% Z(:,9) = isrepeat_same(1:length(Linit.start_time));
% Z(:,10) = isrepeat_different(1:length(Linit.start_time));


Zholder = Z;

% figure;
% scatter(intersect(find(Z(:,2)==1),find(Z(:,1)==1)),ones(1,length(intersect(find(Z(:,2)==1),find(Z(:,1)==1)))),'gsq');
% hold on;
% scatter(intersect(find(Z(:,2)==1),find(Z(:,1)==0)),zeros(1,length(intersect(find(Z(:,2)==1),find(Z(:,1)==0)))),'rsq');
% hold on;
% scatter(intersect(find(Z(:,2)==2),find(Z(:,1)==1)),ones(1,length(intersect(find(Z(:,2)==2),find(Z(:,1)==1)))),'go');
% hold on;
% scatter(intersect(find(Z(:,2)==2),find(Z(:,1)==0)),zeros(1,length(intersect(find(Z(:,2)==2),find(Z(:,1)==0)))),'ro');
% hold on;
% plot([1:size(Z,1)],Z(:,8));
%%  Alignment Section
% close all
pre =5;
post = 3;
clearTrue = 0;
lasShowVal = 1;
baseEnds = [-0.5 0];

filtSize = [5 2];
%Z = Zholder(find(Zholder(:,3) <600),:);

% Define the duration of the delay period
stimLength = 0.6;

visCorrL = intersect(intersect(intersect(intersect(intersect(find(Z(:,1) == 1), find(Z(:,2) == 1)),  find(Z(:,5) == 0)), find(Z(:,8)>=0)), find(Z(:,6)>=0)), [laserinit:length(Z(:,1))]);
audCorrL = intersect(intersect(intersect(intersect(intersect(find(Z(:,1) == 1), find(Z(:,2) == 2)),  find(Z(:,5) == 0)), find(Z(:,8)>=0)), find(Z(:,6)>=0)), [laserinit:length(Z(:,1))]);
visCorrR = intersect(intersect(intersect(intersect(intersect(find(Z(:,1) == 1), find(Z(:,2) == 1)),  find(Z(:,5) == 1)), find(Z(:,8)>=0)), find(Z(:,6)>=0)), [laserinit:length(Z(:,1))]);
audCorrR = intersect(intersect(intersect(intersect(intersect(find(Z(:,1) == 1), find(Z(:,2) == 2)),  find(Z(:,5) == 1)), find(Z(:,8)>=0)), find(Z(:,6)>=0)), [laserinit:length(Z(:,1))]);
% visCorrP = intersect(intersect(find(Z(:,1) == 1), find(Z(:,2) == 1)), find(Z(:,8)>=0.65));
% audCorrP = intersect(intersect(find(Z(:,1) == 1), find(Z(:,2) == 2)), find(Z(:,8)>=0.65));
visCorr = [visCorrR;visCorrL];
audCorr = [audCorrR;audCorrL];
visP=Z(visCorr,8);
audP=Z(audCorr,8);

includeSet = [1:size(num_seq,1)]%14:16]%3 6 17 21 26 28 29];

spikeRangeThresh = 0.03;

for r = 1%:lasShowVal
% visCorr = [visCorrR;visCorrL];
% audCorr = [audCorrR;audCorrL];
% % 
%  visCorr = visCorr(find(visCorr > 40));
%  audCorr = audCorr(find(audCorr > 40));
    if r == 1
        visCorrBase = visCorr;
        if lasShowVal > 1
            visCorr = intersect(visCorr,find(Z(:,6) == 0));
            audCorr = intersect(audCorr,find(Z(:,6) == 0));
        end
        audCorrBase = audCorr;
        spikeRangeThresh = spikeRangeThresh;
        postTag = 'baseline';
    else
        visCorr = intersect(visCorr,find(Z(:,6) == 1));%intersect(find(Z(:,2) == 1),find(Z(:,6) == 1));%
        visCorrLas = visCorr;
        audCorr = intersect(audCorr,find(Z(:,6) == 1));%intersect(find(Z(:,2) == 2),find(Z(:,6) == 1));%
        audCorrLas = audCorr;
        spikeRangeThresh = spikeRangeThresh+0.005;
        postTag = 'laser';
    end
        

mixCorr = [visCorr;audCorr];
sortVecTrials = randperm(length(mixCorr))';

mixCorr = [sortVecTrials mixCorr];
mixCorr = sortrows(mixCorr);
mixCorr = mixCorr(2:round(length(mixCorr)/2),2);



strapLength = 2;

winS = 0;
winEnd = stimLength;

bin = 0.025;
filtWidth = 0.1;

for nn = includeSet;
   handle2=figure(2);
   hold on 
    for m = 1:4
       if m == 1
           f = visCorr;
           stringVal = 'Visual';
           stringVals{m} = stringVal;
       elseif m == 2
           f = audCorr;
           stringVal = 'Auditory';
           stringVals{m} = stringVal;
       elseif m == 3
           f = mixCorr;
           stringVal = 'Mixed';
           stringVals{m} = stringVal;
       else
           f = mixCorr;
           stringVal = 'Baseline';
           stringVals{m} = stringVal;
           winS = -stimLength;
           winEnd = 0; 
       end
        eval(['cl' num2str(nn) '_align = align(cl' num2str(nn) ',Linit,f,pre,post);'])
        eval(['clcurr_align = cl' num2str(nn) '_align;'])
        eval(['clcurr_align' num2str(m) ' = cl' num2str(nn) '_align;']) % save the raster spike timing CG

        if m < 3
%             subplot(2,2,(2*m-1));
            %plot(clcurr_align,'hist','k',30);title([num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2))])
        end
        
        timestamps = clcurr_align.timestamp;
        trial_id = clcurr_align.trial_id';
        trials = unique(trial_id);
        
        clear clCell_align
        clCell_align = []
        qAdd = 1;
        
        if m < 3
        for q = 1:length(f);
            if ~isempty(intersect(q,trials))
                clCell_align{q} = timestamps(find(trial_id == q));
%                 clCell_align{q} = clCell_align{q}(intersect(find(clCell_align{q} > winS),find(clCell_align{q} < winEnd)));
            else
                clCell_align{q} = [];
            end
            fine = f;
        end
        else
            for q = 1:length(f);
            if ~isempty(intersect(q,trials))
                for ii = 1:strapLength
                    clCell_align{q} = timestamps(find(trial_id == q));
%                     clCell_align{q} = clCell_align{q}(intersect(find(clCell_align{q} > winS),find(clCell_align{q} < winEnd)));
                    if m > 3
                        clCell_align{q} = clCell_align{q}+stimLength;
                    end
                    trialJittered = jitterspike(clCell_align{q},0.1);
                     trialJittered(find(trialJittered < 0)) = trialJittered(find(trialJittered < 0))+0.1;
                     trialJittered(find(trialJittered > 0.5)) = trialJittered(find(trialJittered > 0.5))-0.1;
                    clCell_align{q} = sortrows(trialJittered);
                    qAdd = qAdd+1;
                end
            else
                clCell_align{q} = [];
            end
            end
            fine = 1:length(clCell_align);
            fine = fine';
        end
        eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal ' = clCell_align;'])
%         [spikeMean,distances] = spikeDistEstimate(clCell_align,stimLength);% calculate raw peak timing CG
%         [spikePeak,distaste,peakDistances] = spikePeakEstimate(clCell_align,stimLength);
%         spikeMeanCombo{m} = spikeMean;
%         peakDistanceCombo{m} = peakDistances;
%         
%         eval(['spikeMean' stringVal '{nn} = spikeMean;'])% save raw peak timing CG
        
       
%         subplot(2,2,m+4);plot([spikeMean;spikeMean],[ones(1,length(spikeMean));2*ones(1,length(spikeMean))],'k')
%         xlim([-pre post])
        
%         sortVec = [distances' fine];
%         sortVec = sortrows(sortVec);
%         
        
        qfire = [];
        for q = 1:length(f);
            if ~isempty(intersect(q,trials))
                clCell_align{q} = timestamps(find(trial_id == q));
                curr = clCell_align{q};
                if (length(intersect(find(curr < winEnd),find(curr > winS))) > 0)% delete delay null spike trial CG
                    qfire = [qfire q];
                end
            end
        end
        clear clCell_align
        clCell_align = []
        qAdd = 1;
        if m < 3
        for q = qfire;
            qAdd = qAdd+1;
            if ~isempty(intersect(q,trials))
                clCell_align{qAdd} = timestamps(find(trial_id == q));
            end
        end
        else
            for q = qfire;
            if ~isempty(intersect(q,trials))
                for ii = 1
                    trialJittered = jitterspike(timestamps(find(trial_id == q)),0.15);
                    if m > 3
                    trialJittered = trialJittered+stimLength;
                    end
                    trialJittered(find(trialJittered < 0)) = trialJittered(find(trialJittered < 0))+0.15;
                    trialJittered(find(trialJittered > 0.5)) = trialJittered(find(trialJittered > 0.5))-0.15;
                    clCell_align{qAdd} = sortrows(trialJittered);
                    qAdd = qAdd+1;
                end
            else
                clCell_align{q} = [];
            end
            end 
        end
%         if size(clCell_align) >0
%             [pairWise,pairMat] = pairwiseComp(clCell_align,stimLength);
%         
%             basePosibilities(m) = mean(pairWise);
%         
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'pairComp = pairWise;'])
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'pairMat = pairMat;'])
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'Distances = distances;'])
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'SpikeMean = spikeMean;'])
%         else
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'pairComp = 0;'])
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'pairMat = 0;'])
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'Distances = 0;'])
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'SpikeMean = 0;'])
%         end
        
        
        if m < 4
%             hold on;subplot(3,2,(2*m));plot([spikeMean;spikeMean],[zeros(1,length(spikeMean));(length(f)+1)*ones(1,length(spikeMean))],'y')
            %hold on;subplot(3,2,(2*m));plot(clcurr_align,'raster')
        else
%             basePosibilities = basePosibilities(1:2);
%             basePos = find(basePosibilities == min(basePosibilities));
%             hold on;subplot(3,2,(2*basePos));plot([spikeMean-stimLength;spikeMean-stimLength],[zeros(1,length(spikeMean));(length(f)+1)*ones(1,length(spikeMean))],'c')
        end
        eval(['cl' num2str(nn) '_align = align(cl' num2str(nn) ',Linit,f,pre,post);'])
    end
    
    for m = 1:2
       if m == 1
           f = visCorr;
           stringVal = 'Visual';
           stringVals{m} = stringVal;
           ZP=visP;
       elseif m == 2
           f = audCorr;
           stringVal = 'Auditory';
           stringVals{m} = stringVal;
           ZP=audP;
       end
%         clear currSpikeMeanRef
%         currSpikeMean = spikeMeanCombo{m};
%         currPeakDistances = peakDistanceCombo{m};
%         counter = 0;
%         lineSize = length(currPeakDistances{1})+1;
%         for i = 1:length(currSpikeMean)
%             if m < 3
%                 spikeDistThreshComb = min([abs(nanmedian(currPeakDistances{i})) abs(nanmean(currPeakDistances{i}))])
%                 spikeDistThreshTest(nn,i) = spikeDistThreshComb;
%                 eval(['spikeDistThreshTest' num2str(m) '(nn,i) = spikeDistThreshComb;'])%save raw spikeDist CG
%             else
%                 spikeDistThreshComb = 2*abs(nanmedian(currPeakDistances{i}))
%             end
%             if spikeDistThreshComb < spikeRangeThresh % find spike dist< thresh peak timing red line CG
%                 currSpikeMeanRef(i) = currSpikeMean(i)
%                 counter = counter+1;
%             end
%         end
%         hold on;subplot(2,2,(2*m));plot([currSpikeMean;currSpikeMean],[zeros(1,length(currSpikeMean));(lineSize)*ones(1,length(currSpikeMean))],'y')
%         if counter > 0
%             if counter < 15
%         	currSpikeMeanRef = currSpikeMeanRef(find(currSpikeMeanRef > 0));
%             hold on;subplot(2,2,(2*m));plot([currSpikeMeanRef;currSpikeMeanRef],[zeros(1,length(currSpikeMeanRef));(lineSize)*ones(1,length(currSpikeMeanRef))],'r')
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'SpikePeak = currSpikeMeanRef;']) % save peak timing CG
%             end
%         else
%             eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal 'SpikePeak = NaN;'])
%         end
        
        eval(['clcurr_align = clcurr_align' num2str(m) ';'])
        %
%         if exist('currSpikeMeanRef')
%         closeSpikesL = intersect(find(clcurr_align.timestamp > (max(currSpikeMeanRef)-0.065)),find(clcurr_align.timestamp < (max(currSpikeMeanRef))))
%         closeSpikesH = intersect(find(clcurr_align.timestamp < (max(currSpikeMeanRef)+0.065)),find(clcurr_align.timestamp > (max(currSpikeMeanRef))))
%         clcurr_align.timestamp(closeSpikesL) = clcurr_align.timestamp(closeSpikesL)+0.045;
%         clcurr_align.timestamp(closeSpikesH) = clcurr_align.timestamp(closeSpikesH)-0.045;
%         end
        %
%         hold on;subplot(2,2,(2*m));plot(clcurr_align,'raster')
%         ax1 = gca;
%         ax1_pos = ax1.Position; % position of first axes
%         ax2 = axes('Position',ax1_pos,...
%             'XAxisLocation','top',...
%             'YAxisLocation','right',...
%             'Color','none');
%         plot(ZP,[1:length(ZP)],'Parent',ax2,'Color','r')
        timestamps = clcurr_align.timestamp;
        trial_id = clcurr_align.trial_id';
        trials = unique(trial_id);
        
        clear clCell_align
        
        qAdd = 1;
        
        if m < 3
        for q = 1:length(f);
            if ~isempty(intersect(q,trials))
                clCell_align{q} = timestamps(find(trial_id == q));
%                 clCell_align{q} = clCell_align{q}(intersect(find(clCell_align{q} > winS),find(clCell_align{q} < winEnd)));
            else
                clCell_align{q} = [];
            end
            fine = f;
        end
         if m < 3
%             [spikeRate,error,errorBounds,spikeMat,time] = spikeRateEst(clCell_align,bin,filtWidth,pre,post);
%             subplot(2,2,(2*m-1));
%             spikePlot = nanconv(spikeRate,gausswin(5)'./sum(gausswin(5)))
%             plot(time,spikePlot,'k');title(['Unit' num2str(nn) 'ID' num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2))]);hold on
%             errorBoundsPlot = [spikePlot+mean(error)/1.5;spikePlot-mean(error)/1.5]
%             plot(time,errorBoundsPlot(1,:),'Color',[0.5 0.5 0.5]); plot(time,errorBoundsPlot(2,:),'Color',[0.5 0.5 0.5])
         end
      
            %[spikeRate,error,errorBounds,spikeMat,time] = spikeRateEst(clCell_align,bin,filtWidth,pre,post);
            [zscore,firingRate,baseMean,peakScore,time] = spikeZscore(clCell_align,bin,baseEnds,pre,post,0.1,[]);
%             subplot(2,2,(2*m-1));
            zscore = nanconv(0.25+zscore,gausswin(filtSize(1),filtSize(2))'./sum(gausswin(filtSize(1),filtSize(2))));
%             hold off
%             plot(time,1.5*zscore,'k');title(['Unit' num2str(nn) 'ID' num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2))]);hold on 
%             errorBoundsAlt = 
%             plot(time,errorBoundsAlt(:,1),'Color',[0.5 0.5 0.5]); plot(time,errorBoundsAlt(:,2),'Color',[0.5 0.5 0.5])
        end
    end
    if r == 1
%         cd(figpath);
%         saveas(gcf,['Example AMP4 S' num2str(session_num) ' U' num2str(nn) 'ID' num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2)) ' Tuning Combined.fig'])
%         saveas(gcf,['Example AMP4 S' num2str(session_num) ' U' num2str(nn) 'ID' num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2)) ' Tuning Combined.png'])
%         close all;
    else
%         saveas(gcf,['Example WT19 S' num2str(session_num) ' U' num2str(nn) 'ID' num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2)) '  MDSSFO Tuning.fig'])
%         close 
    end
    if clearTrue == 1
    clear(['cl' num2str(nn) ''],['cl' num2str(nn) '_align']) 
    end
    %eval(['[h,pVal] =kstest2(cl' num2str(nn) '_cellData' postTag '.VisualDistances,cl' num2str(nn) '_cellData' postTag '.MixedDistances);'])
    %pTest(n) = pVal;
end
%pTest = pTest';

clear Sc Se clcurr_aligncl clcurr_align clcurr_align1 clcurr_align2 clcurr_align3 clcurr_align4
end

%  Alignment Section  Incorrect Trials

close all
pre = 1.5;
post = 1.5;
clearTrue = 0;
lasShowVal = 2;

visIncL = intersect(intersect(find(Z(:,1) == 0), find(Z(:,2) == 1)),  find(Z(:,5) == 0));
audIncL = intersect(intersect(find(Z(:,1) == 0), find(Z(:,2) == 2)),  find(Z(:,5) == 0));
visIncR = intersect(intersect(find(Z(:,1) == 0), find(Z(:,2) == 1)),  find(Z(:,5) == 1));
audIncR = intersect(intersect(find(Z(:,1) == 0), find(Z(:,2) == 2)),  find(Z(:,5) == 1));
%includeSet = [2:29];


spikeRangeThresh = 0.025;

for r = 1%:lasShowVal
visInc = [visIncR;visIncL];
audInc = [audIncR;audIncL];
    if r == 1
        %visInc = intersect(visInc,find(Z(:,6) == 0));
        %audInc = intersect(audInc,find(Z(:,6) == 0));
        spikeRangeThresh = spikeRangeThresh;
        postTag = 'baseline';
    else
        visInc = intersect(visInc,find(Z(:,6) == 1));
        audInc = intersect(audInc,find(Z(:,6) == 1));
        spikeRangeThresh = spikeRangeThresh-0.005;
        postTag = 'laser';
        %clearTrue = 1;
    end
        

mixInc = [visInc;audInc];
sortVecTrials = randperm(length(mixInc))';

mixInc = [sortVecTrials mixInc];
mixInc = sortrows(mixInc);
mixInc = mixInc(2:round(length(mixInc)/2),2);

stimLength = 0.6;

strapLength = 2;

winS = 0;
winEnd = stimLength;

bin = 0.025;
filtWidth = 0.18;

for nn = includeSet;
%    figure
%    hold on 
    for m = 1:4
       if m == 1
           f = visInc;
           stringVal = 'Visual';
           stringVals{m} = stringVal;
       elseif m == 2
           f = audInc;
           stringVal = 'Auditory';
           stringVals{m} = stringVal;
       elseif m == 3
           f = mixInc;
           stringVal = 'Mixed';
           stringVals{m} = stringVal;
       else
           f = mixInc;
           stringVal = 'Baseline';
           stringVals{m} = stringVal;
           winS = -stimLength;
           winEnd = 0; 
       end
        eval(['cl' num2str(nn) '_alignInc = align(cl' num2str(nn) ',Linit,f,pre,post);'])
        eval(['clcurr_align = cl' num2str(nn) '_alignInc;'])
        eval(['clcurr_align' num2str(m) 'Inc = cl' num2str(nn) '_alignInc;'])

        if m < 3
%             subplot(2,2,(2*m-1));
            %plot(clcurr_align,'hist','k',30);title([num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2))])
        end
        
        timestamps = clcurr_align.timestamp;
        trial_id = clcurr_align.trial_id';
        trials = unique(trial_id);
        
        clear clCell_align
        clCell_align = []
        qAdd = 1;
        
        if m < 3
        for q = 1:length(f);
            if ~isempty(intersect(q,trials))
                clCell_align{q} = timestamps(find(trial_id == q));
%                 clCell_align{q} = clCell_align{q}(intersect(find(clCell_align{q} > winS),find(clCell_align{q} < winEnd)));
            else
                clCell_align{q} = [];
            end
            fine = f;
        end
        else
            for q = 1:length(f);
            if ~isempty(intersect(q,trials))
                for ii = 1:strapLength
                    clCell_align{q} = timestamps(find(trial_id == q));
%                     clCell_align{q} = clCell_align{q}(intersect(find(clCell_align{q} > winS),find(clCell_align{q} < winEnd)));
                    if m > 3
                        clCell_align{q} = clCell_align{q}+stimLength;
                    end
                    trialJittered = jitterspike(clCell_align{q},0.1);
                     trialJittered(find(trialJittered < 0)) = trialJittered(find(trialJittered < 0))+0.1;
                     trialJittered(find(trialJittered > 0.5)) = trialJittered(find(trialJittered > 0.5))-0.1;
                    clCell_align{q} = sortrows(trialJittered);
                    qAdd = qAdd+1;
                end
            else
                clCell_align{q} = [];
            end
            end
            fine = 1:length(clCell_align);
            fine = fine';
        end
        eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal ' = clCell_align;'])
%         [spikeMean,distances] = spikeDistEstimate(clCell_align,stimLength);
%         [spikePeak,distaste,peakDistances] = spikePeakEstimate(clCell_align,stimLength);
%         spikeMeanCombo{m} = spikeMean;
%         peakDistanceCombo{m} = peakDistances;
%         
%         eval(['spikeMean' stringVal '{nn} = spikeMean;'])
        
        if m < 3
%             [spikeRate,error,errorBounds,spikeMat,time] = spikeRateEst(clCell_align,bin,filtWidth,pre,post);
%             subplot(2,2,(2*m-1));
%             spikePlot = nanconv(spikeRate,[0.15 0.25 1 0.25 0.15]/2);
% %             plot(time,spikePlot,'k');title([num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2))]);hold on
%             errorBoundsPlot = [spikePlot+mean(error);spikePlot-mean(error)];
% %             plot(time,errorBoundsPlot(1,:),'Color',[0.5 0.5 0.5]); plot(time,errorBoundsPlot(2,:),'Color',[0.5 0.5 0.5])
        end
%         subplot(2,2,m+4);plot([spikeMean;spikeMean],[ones(1,length(spikeMean));2*ones(1,length(spikeMean))],'k')
%         xlim([-pre post])
        
%         sortVec = [distances' fine];
%         sortVec = sortrows(sortVec);
%         
        
        qfire = [];
        for q = 1:length(f);
            if ~isempty(intersect(q,trials))
                clCell_align{q} = timestamps(find(trial_id == q));
                curr = clCell_align{q};
                if (length(intersect(find(curr < winEnd),find(curr > winS))) > 0)
                    qfire = [qfire q];
                end
            end
        end
        clear clCell_align
        clCell_align = []
        qAdd = 1;
        if m < 3
        for q = qfire;
            qAdd = qAdd+1;
            if ~isempty(intersect(q,trials))
                clCell_align{qAdd} = timestamps(find(trial_id == q));
            end
        end
        else
            for q = qfire;
            if ~isempty(intersect(q,trials))
                for ii = 1
                    trialJittered = jitterspike(timestamps(find(trial_id == q)),0.15);
                    if m > 3
                    trialJittered = trialJittered+stimLength;
                    end
                    trialJittered(find(trialJittered < 0)) = trialJittered(find(trialJittered < 0))+0.15;
                    trialJittered(find(trialJittered > 0.5)) = trialJittered(find(trialJittered > 0.5))-0.15;
                    clCell_align{qAdd} = sortrows(trialJittered);
                    qAdd = qAdd+1;
                end
            else
                clCell_align{q} = [];
            end
            end 
        end
%         if size(clCell_align) >0
%             [pairWise,pairMat] = pairwiseComp(clCell_align,stimLength);
%         
%             basePosibilities(m) = mean(pairWise);
%         
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'pairComp = pairWise;'])
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'pairMat = pairMat;'])
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'Distances = distances;'])
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'SpikeMean = spikeMean;'])
%         else
%             basePosibilities=[];
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'pairComp = 0;'])
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'pairMat = 0;'])
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'Distances = 0;'])
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'SpikeMean = 0;'])
%         end
        
        
%         if m < 4
% %             hold on;subplot(3,2,(2*m));plot([spikeMean;spikeMean],[zeros(1,length(spikeMean));(length(f)+1)*ones(1,length(spikeMean))],'y')
%             %hold on;subplot(3,2,(2*m));plot(clcurr_align,'raster')
%         else
%             if length(basePosibilities)>0
%             basePosibilities = basePosibilities(1:2);
%             basePos = find(basePosibilities == min(basePosibilities));
% %             hold on;subplot(3,2,(2*basePos));plot([spikeMean-stimLength;spikeMean-stimLength],[zeros(1,length(spikeMean));(length(f)+1)*ones(1,length(spikeMean))],'c')
%             end
%         end
        eval(['cl' num2str(nn) '_align = align(cl' num2str(nn) ',Linit,f,pre,post);'])
    end
    
    for m = 1:2
       if m == 1
           f = visInc;
           stringVal = 'Visual';
           stringVals{m} = stringVal;
       elseif m == 2
           f = audInc;
           stringVal = 'Auditory';
           stringVals{m} = stringVal;
       end
%         clear currSpikeMeanRef
%         currSpikeMean = spikeMeanCombo{m};
%         currPeakDistances = peakDistanceCombo{m};
%         counter = 0;
%         lineSize = length(currPeakDistances{1})+1;
%         for i = 1:length(currSpikeMean)
%             if m < 3
%                 spikeDistThreshComb = min([abs(nanmedian(currPeakDistances{i})) abs(nanmean(currPeakDistances{i}))])
%                 spikeDistThreshTest(nn,i) = spikeDistThreshComb;
%                 eval(['spikeDistThreshTestInc' num2str(m) '(nn,i) = spikeDistThreshComb;'])
%             else
%                 spikeDistThreshComb = 2*abs(nanmedian(currPeakDistances{i}))
%             end
%             if spikeDistThreshComb < spikeRangeThresh
%                 currSpikeMeanRef(i) = currSpikeMean(i)
%                 counter = counter+1;
%             end
%         end
%         hold on;subplot(2,2,(2*m));plot([currSpikeMean;currSpikeMean],[zeros(1,length(currSpikeMean));(lineSize)*ones(1,length(currSpikeMean))],'y')
%         if counter > 0
%             if counter < 15
%         	currSpikeMeanRef = currSpikeMeanRef(find(currSpikeMeanRef > 0));
% %             hold on;subplot(2,2,(2*m));plot([currSpikeMeanRef;currSpikeMeanRef],[zeros(1,length(currSpikeMeanRef));(lineSize)*ones(1,length(currSpikeMeanRef))],'r')
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'SpikePeak = currSpikeMeanRef;'])
%             end
%         else
%             eval(['cl' num2str(nn) '_cellData' postTag 'Inc.' stringVal 'SpikePeak = NaN;'])
%         end
        eval(['clcurr_align = clcurr_align' num2str(m) 'Inc;'])
%         hold on;subplot(2,2,(2*m));plot(clcurr_align,'raster')
    end    
%     eval(['clcurr_align = clcurr_align3;'])
%     hold on;subplot(3,2,(6));plot(clcurr_align,'raster')
    
%     figure
%     eval(['cdfplot(cl' num2str(nn) '_cellData' postTag '.VisualpairComp);']); hold on 
%     eval(['cdfplot(cl' num2str(nn) '_cellData' postTag '.AuditorypairComp);'])
%     eval(['cdfplot(cl' num2str(nn) '_cellData' postTag '.MixedpairComp);'])
%     eval(['cdfplot(cl' num2str(nn) '_cellData' postTag '.BaselinepairComp);'])
%     title([num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2))])
    
%     figure
%     eval(['cdfplot(cl' num2str(nn) '_cellData' postTag '.VisualDistances);']); hold on 
%     eval(['cdfplot(cl' num2str(nn) '_cellData' postTag '.AuditoryDistances);'])
%     eval(['cdfplot(cl' num2str(nn) '_cellData' postTag '.MixedDistances);'])
%     title([num2str(num_seq(nn,1)),'-',num2str(num_seq(nn,2))])
    if clearTrue == 1
    clear(['cl' num2str(nn) ''],['cl' num2str(nn) '_alignInc']) 
    end
    %eval(['[hInc,pValInc] =kstest2(cl' num2str(nn) '_cellData' postTag '.VisualDistances,cl' num2str(nn) '_cellData' postTag '.MixedDistances);'])
    %pTest(n) = pVal;
end
%pTest = pTest';
   
clear Sc Se clcurr_aligncl clcurr_align clcurr_align1 clcurr_align2 clcurr_align3 clcurr_align4 
end
cd(figpath);

if extt == 0
save(strcat(Fpath(end-18:end),'_PFCMDSection',num2str(cn),'TT',num2str(tn),'.mat'),'-v7.3');
elseif extt == 1
    save(strcat(Fpath(4:end),'_PFCMDSection',num2str(cn),'TT',num2str(tn),'.mat'),'-v7.3');
end
% save(strcat(Fpath(end-18:end),'_PFCMDSection',num2str(cn),'P.mat'),'visCorr','audCorr','sn');
clearvars -except sn cn tn
close all;
        end
    end
end
%%
% %Unit creator
% for n = 1:(length(num_seq(:,1)))
% eval(['unit' num2str(n) ' = cl' num2str(n) '.timestamp;'])
% eval(['unit' num2str(n) ' = unit' num2str(n) ';'])
% eval(['clear cl' num2str(n) ';'])
% eval(['clear cl' num2str(n) '_align;'])
% eval(['clear cl' num2str(n) '_alignInc;'])
% end
%%
if extt == 0
    save(strcat(Fpath(end-18:end),'_PFCMDSection',num2str(cn),'.mat'),'-v7.3');
elseif extt ==1
    save(strcat(Fpath,'_PFCMDSection',num2str(cn),'.mat'),'-v7.3');
end;
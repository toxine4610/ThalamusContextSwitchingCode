
function BuildData_Shrew


%%

%%%==== process accumulation of evidence experiments
clear all;
MouseID = 'Peach';
DataLocation = '/Volumes/Data/ShrewPeach/2018-01-29/';
behaviorfilename.C1 = [ DataLocation filesep 'Zelda_20180129.txt'];
behaviorfilename.C2 = behaviorfilename.C1;

cd(DataLocation);
addpath(DataLocation);

switch MouseID
    case 'SOMCRE1'
        TT_to_process = 14;
    case 'WT24'
        TT_to_process = [1:15,17:24];
    case 'Constantine'
        TT_to_process = [1:15,18];
    case 'BPM1_3';
        TT_to_process =  [1:24];
    case 'Accum'
        TT_to_process = [1:16];
    case 'Peach'
        TT_to_process = 10:32 % setdiff( [1:32],[6,7,8] );
end;


% addpath('G:\Session Import');
% addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\ForShrew_curveFit')
% addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\Rajeev_Code');
% addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\Rajeev_Code\Clustering and Basic Analysis');
% addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\Rajeev_Code\Clustering and Basic Analysis\mclust-3.4');
addpath('/Volumes/Data/Session Import');
addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForShrew_curveFit')
addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/Rajeev_Code/');
addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/Rajeev_Code/Clustering and Basic Analysis/');
addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/Rajeev_Code/Clustering and Basic Analysis/mclust-3.4/');


%%
for tn = TT_to_process
    
    %     if getappdata(h,'canceling')
    %         break
    %     end
    
    %% Make a Temp Session
    extt = 1; % save to drobo;
    
    session_num =  1;
    
    switch MouseID
        case 'SOMCRE1'
            Se = sessions();
            Se = add(Se,'CMO',...
                [2014 04 13 15 11 22],...
                DataLocation,...
                'autoCluster',32,...
                'File 6',...
                [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4],... % channels per electrode/home
                [NaN NaN NaN NaN NaN NaN]);
        case 'WT24'
            Se = sessions();
            Se = add(Se,'CMO',...
                [2014 04 13 15 11 22],...
                DataLocation,...
                'autoCluster',32,...
                'File 6',...
                [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 ],... % channels per electrode/home
                [NaN NaN NaN NaN NaN NaN]);
        case 'BPM1_3'
            Se = sessions();
            Se = add(Se,'CMO',...
                [2014 04 13 15 11 22],...
                DataLocation,...
                'autoCluster',32,...
                'File 6',...
                [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 ],... % channels per electrode/home
                [NaN NaN NaN NaN NaN NaN]);
        case 'Constantine'
            Se = sessions();
            Se = add(Se,'CMO',...
                [2014 04 13 15 11 22],...
                DataLocation,...
                'autoCluster',32,...
                'File 6',...
                [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4],... % channels per electrode/home
                [NaN NaN NaN NaN NaN NaN]);
            
        case 'Accum'
            Se = sessions();
            Se = add(Se,'CMO',...
                [2014 04 13 15 11 22],...
                DataLocation,...
                'autoCluster',32,...
                'File 6',...
                [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4],... % channels per electrode/home
                [NaN NaN NaN NaN NaN NaN]);
            
        case 'Peach'
            Se = sessions();
            Se = add(Se,'CMO',...
                [2014 04 13 15 11 22],...
                DataLocation,...
                'Neuralynx G',32,...
                'File 6',...
                [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4],... % channels per electrode/home
                [NaN NaN NaN NaN NaN NaN]);
    end
    
    Fpath = Se.folder{session_num};
    figpath=strcat(Fpath, filesep, 'GlobalpeakC1');
    mkdir(figpath);
    
    %% Parse behavior text file.
    
    %========== Load TTL Values
    eventData = read_cheetah_data([Se.folder{session_num} filesep 'Events.nev']);
    TTLvalues = unique( eventData.TTLval );
    switch MouseID
        case 'SOMCRE1'
            TTLbehave =  1;
            TTstim1   =  4;
        case 'WT24'
            TTLbehave =  128;
            TTLstim1  =  64;
        case 'Constantine'
            TTLbehave = 2;
            TTLstim1  =  4;
        case 'BPM1_3'
            TTLbehave = 2;
            TTLstim1  =  4;
        case 'Accum'
            TTLbehave =  1;
            TTstim1   =  4;
            
        case 'Peach'
            TTLbehave = 2;
            TTstim1   =  4;
    end;
    
    TTLbehave2 = 193;
    
    
    eventData = read_cheetah_data([Se.folder{session_num}  filesep  'Events.nev']);
    eventDataholder = eventData;
    % Lstim.start_time = eventData.ts(find(eventData.TTLval == TTLstim1))';
    % Lstim.pulse_width = diff(eventData.ts);
    % Lstim.pulse_width = Lstim.pulse_width(find(eventData.TTLval == TTLstim1))';
    L.start_time = eventData.ts(find((eventData.TTLval == TTLbehave)+(eventData.TTLval == TTLbehave2)))';
    
    
    %%
    correctFlag = 1;
    
    Out = parseShrewBehavior(behaviorfilename.C1);
    D = Out.D;
    Dcorrected = Out.Dcorrected;
    
    %%
    
    initT = Out.initT;
    % if Pos offset = -2 : start from 3rd timestamp (missed 2)
        initT = Out.initT(1:end-19);
        Out.ts  = Out.ts(1:end-19);
    
    
    
    [a,b,offset] = alignsignals(round(diff(Out.ts/1000)),round(diff(L.start_time)));
    [a,b,negOff] = alignsignals(round(diff(wrev(Out.ts)/1000)),round(diff(wrev(L.start_time))));
    
    disp(['Pos Offset = ' num2str(offset) '  Neg Offset = ' num2str(negOff)]);
    
    
    %%
    if negOff >= 0 && offset <= 0
        
        testDiff = L.start_time(1:(length(L.start_time)-negOff)) - initT((-offset+1):(length(initT)))'/1000;
        figure(1);
        subplot(2,1,1);plot(diff(testDiff))
        subplot(2,1,2);plot(diff(Out.ts))
        
        drawnow;
        Linit.start_time = L.start_time(1:(length(L.start_time)-negOff))  - initT((-offset+1):(length(initT)))'/1000;
        
        if correctFlag == 1
            Linit.start_time = Linit.start_time( Out.indexOK);
            D = Out.Dcorrected;
        end;
        
    end;
    
    
    if negOff < 0 && offset == 0
        Linit.start_time = L.start_time((offset+1):length(L.start_time)) - initT((1):(length(initT)+negOff))'/1000;
        testDiff = L.start_time((offset+1):length(L.start_time)) - initT(1:(length(initT)+negOff))'/1000;
        
        figure(1);
        subplot(2,1,1);plot(diff(testDiff))
        subplot(2,1,2);plot(diff(Out.ts))
        
        
%         indexOK = Out.indexOK(1:end-1);
%         
%         Linit.start_time = Linit.start_time(indexOK);
%         D = Out.D;
%         D = D(indexOK,:);
    end;
    
    
    %% Get data from Neuralynx
    autoclust =1
    
    
    
    if autoclust == 1
    
    unitCount = 1;
    
    eventData = read_cheetah_data([Se.folder{session_num} filesep 'Events.nev']);
    
    autoClusterSpikeFiles = dir([ DataLocation filesep 'TT*cluster.mat' ] );
    
    for n = 1:length(autoClusterSpikeFiles)
        TTname = autoClusterSpikeFiles(n).name;
        numTT(n) = str2num(TTname(3:(length(TTname)-11)));
    end
    
    [TTNums,TTOrder] = sort(numTT);
    
    
    %     waitbar(tn/length(TT_to_process),hW, sprintf( ['Processing -- Tet Number: ' '%d', num2str( TTnums( find(TTNums == tn) ) )]  ) );
    
    for i = tn % 1:length(TTNums)
        ind = find( TTNums == tn );
        fname = [ DataLocation filesep 'TT' num2str(TTNums(ind)) '.ntt']
        tetrodeSpikeTimes = read_cheetah_data(fname);
        load([autoClusterSpikeFiles(TTOrder(ind)).name])
        eval(['tetrodeSet = TT' num2str(TTNums(ind)) '.labels;'])
        
        
        Sc = spikesAuto(Se,session_num,TTNums(ind));
        for ii = 1:length(tetrodeSet)
            if max(tetrodeSet{ii}) < length(tetrodeSpikeTimes.ts)
                cl_holder = clusterAuto(Sc,ii);
                unitSpikeTimes = tetrodeSpikeTimes.ts(tetrodeSet{ii});
                eval(['unitAC' num2str(unitCount) ' = unitSpikeTimes;'])
                eval(['cl' num2str(unitCount) ' = cl_holder;'])
                autoClusterSpikeTimes{i}{ii} = unitSpikeTimes;
                num_seq(unitCount,:) = [(TTNums(ind)) ii];
            else
                eval(['tetrodeSpikeTimesAuto = TT' num2str(TTNums(ind)) '.timestamps;'])
                cl_holder = clusterAuto(Sc,ii);
                unitSpikeTimes = tetrodeSpikeTimesAuto(tetrodeSet{ii});
                eval(['unitAC' num2str(unitCount) ' = unitSpikeTimes;'])
                eval(['cl' num2str(unitCount) ' = cl_holder;'])
                autoClusterSpikeTimes{i}{ii} = unitSpikeTimes;
                num_seq(unitCount,:) = [(TTNums(ind)) ii];
            end
            unitCount = unitCount+1;
        end
    end
    
    end;
    
    
%     %% Manual Cluster 
%     
%     if autoclust == 0
%         
%         names = dir('*.cut');
%         all_tt_nums = [];
%         
%         
%         for n = 1:length(names)
%             nameSU = names(n).name;
%             nameSU(strfind(nameSU,'.cut'):length(nameSU)) = [];
%             nameSU(1:2) = [];
%             all_tt_nums = [all_tt_nums str2num(nameSU)];
%         end
%         
%         all_tt_nums = sort(all_tt_nums); %Automatic selection of electrode sets to plot
%         
%         
%         %File Name Adapter
%         
%         names = dir('*.ntt')
%         newNames = names;
%         
%         for n = 1:length(names)
%             newNames(n).name(1) = 'T';
%             newNames(n).name(2) = 'T';
%             dos(['rename "' names(n).name '" "' newNames(n).name '"']);
%         end
%         
%         names = dir('*.cut')
%         newNames = names;
%         
%         for n = 1:length(names)
%             newNames(n).name(1) = 'T';
%             newNames(n).name(2) = 'T';
%             dos(['rename "' names(n).name '" "' newNames(n).name '"']);
%         end
%         
%         all_tt_nums = [];
%         
%         % == manual import
%         fprintf('Extracting Data....');
%         
%         i = 0;
%         
%         for n = 1:(length(all_tt_nums))
%             curr_tt_num = all_tt_nums(n);
%             Sc = spikes(Se,session_num,curr_tt_num);
%             is_cluster = 1;
%             m = 0;
%             while is_cluster == 1
%                 m = m+1;
%                 cl_holder = cluster(Sc,m);
%                 is_full = max(size(cl_holder.timestamp));
%                 if is_full > 1
%                     i = i+1;
%                     eval(sprintf('cl%d = cluster(Sc,m)', i));
%                     num_seq(i,1:2) = [curr_tt_num m];
%                     
%                     
%                     
%                 else
%                     m = m-1;
%                     is_cluster = 0;
%                 end
%             end
%             Sc_unit_count(n,1) = curr_tt_num;
%             Sc_unit_count(n,2) = m;
%         end
%         
%         clear cl_holder i curr_tt_num is_cluster is_full m n
%         fprintf('....Done!\n');
%         
%     end
    
    %% Align spikes
    
    if max(size(tetrodeSet)) == 0;
        disp(['No Units on this Tetrode']);
    end
    
    if max(size(tetrodeSet)) ~= 0
        
        % close all
        pre = 5;
        post = 3;
        clearTrue = 0;
        lasShowVal = 1;
        baseEnds = [-0.5 0];
        filtSize = [5 2];
        postTag = 'baseline';
        
        includeSet = [1:size(num_seq,1)];
        
        % pure correct
        VisL_corr = find( D(:,1) == 1 & D(:,5) == 1 );
        VisR_corr = find( D(:,1) == 1 &  D(:,5) == -1 );
        Vis_corr  = sort( [VisL_corr; VisR_corr] );
        
        % pure correct
        AudL_corr = find( D(:,1) == 1  & D(:,5) == 2 );
        AudR_corr = find( D(:,1) == 1  & D(:,5) == -2 );
        Aud_corr  = sort( [AudL_corr; AudR_corr] );
        
        
        % pure incorrect
        VisL_incorr = find( D(:,1) == 0  & D(:,5) == 1 );
        VisR_incorr = find( D(:,1) == 0  & D(:,5) == -1 );
        Vis_incorr  = sort( [VisL_incorr; VisR_incorr] );
        
        % pure incorrect
        AudL_incorr = find( D(:,1) == 0  & D(:,5) == 2 );
        AudR_incorr = find( D(:,1) == 0  & D(:,5) == -2 );
        Aud_incorr  = sort( [AudL_incorr; AudR_incorr] );
        
        
        for nn = includeSet;
            for m = 1:4
                
                if m == 1
                    f = Vis_corr;
                    stringVal = 'Vis_corr';
                    stringVals{m} = stringVal;
                elseif m == 2
                    f = Aud_corr;
                    stringVal = 'Aud_corr';
                    stringVals{m} = stringVal;
                    
                    
                elseif m == 3
                    f = Vis_incorr;
                    stringVal = 'Vis_incorr';
                    stringVals{m} = stringVal;
                elseif m == 4
                    f = Aud_incorr;
                    stringVal = 'Aud_incorr';
                    stringVals{m} = stringVal;
                    
                end
                eval(['cl' num2str(nn) '_align = align(cl' num2str(nn) ',Linit,f,pre,post);']);
                eval(['clcurr_align = cl' num2str(nn) '_align;'])
                eval(['clcurr_align' num2str(m) ' = cl' num2str(nn) '_align;']) % save the raster spike timing CG
                
                timestamps = clcurr_align.timestamp;
                trial_id = clcurr_align.trial_id';
                trials = unique(trial_id);
                
                clear clCell_align
                clCell_align = []
                qAdd = 1;
                
                for q = 1:length(f);
                    if ~isempty(intersect(q,trials))
                        clCell_align{q} = timestamps(find(trial_id == q));
                        %                 clCell_align{q} = clCell_align{q}(intersect(find(clCell_align{q} > winS),find(clCell_align{q} < winEnd)));
                    else
                        clCell_align{q} = [];
                    end
                    fine = f;
                end
                fine = 1:length(clCell_align);
                fine = fine';
                eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal ' = clCell_align;'])
            end
            
            clear Sc Se clcurr_aligncl clcurr_align clcurr_align1 clcurr_align2 clcurr_align3 clcurr_align4 h
            
        end;
        cd(figpath);
        
        if extt == 0
            save(strcat(Fpath(end-18:end),'_PFCMDSection',num2str(1),'TT',num2str(tn),'.mat'),'-v7.3');
        elseif extt == 1
            save(strcat(date,'_PFCMDSection',num2str(1),'TT',num2str(tn),'.mat'),'-v7.3');
        end
        % save(strcat(Fpath(end-18:end),'_PFCMDSection',num2str(cn),'P.mat'),'AudCorr','audCorr','sn');
        clearvars -except sn cn tn numTT MouseID DataLocation behaviorfilename  TT_to_process h
        close all;
        
    end;
end;

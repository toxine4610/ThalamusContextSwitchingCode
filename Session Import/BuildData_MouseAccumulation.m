
function BuildData_MouseAccumulation


%%
addpath('C:\Users\Halassalab-CG\Desktop\Session Import')
addpath('C:\Users\Halassalab-CG\Desktop\ShrewAnalysis')
%%%==== process accumulation of evidence experiments
clear all;
MouseID = 'BPM4_3';
DataLocation = 'X:\Rajeev\Cheetah\M4-3\2018-08-07_15-55-34';
behaviorfilename.C1 = [ DataLocation filesep 'bpm43_mdterminal_6040_20180807.txt'];
behaviorfilename.C2 = behaviorfilename.C1;

experimentType = 'New';

cd(DataLocation);
addpath(DataLocation);


ValidNum = [];
close
X = dir('*.mat');
for i = 1:size(X,1)
    name = sprintf( X(i).name );
    target1 = strfind(name, 'T');
    target2 = strfind(name, 'c');
    num = name( target1+2 : target2-1  );
    ValidNum = [ValidNum, str2num(num)];
end;

TT_to_process = sort(ValidNum);


if ismac ==1
    addpath('/Volumes/Data/Session Import');
    addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForShrew_curveFit')
    addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/Rajeev_Code/');
    addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/Rajeev_Code/Clustering and Basic Analysis/');
    addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/Rajeev_Code/Clustering and Basic Analysis/mclust-3.4/');
elseif ismac ==0
    addpath('G:\Session Import');
    addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\ForShrew_curveFit');
    addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\Rajeev_Code')
    addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\Rajeev_Code\Clustering and Basic Analysis');
    addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\Rajeev_Code\Clustering and Basic Analysis\mclust-3.4');
    addpath('C:\Users\Halassalab-CG\Dropbox\Rajeev\Rajeev_Code\Clustering and Basic Analysis\readcheetahdata');
end;

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
            
        case 'BPM3_3'
            Se = sessions();
            Se = add(Se,'CMO',...
                [2014 04 13 15 11 22],...
                DataLocation,...
                'autoCluster',32,...
                'File 6',...
                [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 ],... % channels per electrode/home
                [NaN NaN NaN NaN NaN NaN]);
            
        case 'BPM4_3'
            Se = sessions();
            Se = add(Se,'CMO',...
                [2014 04 13 15 11 22],...
                DataLocation,...
                'autoCluster',32,...
                'File 6',...
                [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4],... % channels per electrode/home
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
        case 'BPM3_3'
            TTLbehave = 128;
            TTLstim1  =  64;
        case 'BPM4_3'
            TTLbehave =  1;
            TTstim1   =  0;
    end;
    
    TTLbehave2 = 193;
    
    
    eventData = read_cheetah_data([Se.folder{session_num}  filesep  'Events.nev']);
    eventDataholder = eventData;
    % Lstim.start_time = eventData.ts(find(eventData.TTLval == TTLstim1))';
    % Lstim.pulse_width = diff(eventData.ts);
    % Lstim.pulse_width = Lstim.pulse_width(find(eventData.TTLval == TTLstim1))';
    L.start_time = eventData.ts(find((eventData.TTLval == TTLbehave)+(eventData.TTLval == TTLbehave2)))';
    %     L.start_time = L.start_time(41:end);
    
    %%
    correctFlag = 1;
    
    Out = parse4AFC_shrews_lasers(behaviorfilename.C1,[] );
    D = Out.D;
    Dcorrected = Out.Dcorrected;
    
    %%
    
    initT = 0*Out.initT;
    % if Pos offset = -2 : start from 3rd timestamp (missed 2)
    %         initT = Out.initT(2:end);
    %         Out.ts  = Out.ts(2:end);
    
    
    
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
        
        %         if correctFlag == 1
        %             Linit.start_time = Linit.start_time( Out.indexOK);
        %             D = Out.Dcorrected;
        %         end;
        
    end;
    
    if negOff >= 0 && offset > 0
         Linit.start_time = L.start_time((offset+1):length(L.start_time)) - initT((1):(length(initT)+negOff))'/1000
         testDiff = L.start_time((offset+1):length(L.start_time)) - initT(1:(length(initT)+negOff))'/1000
             figure(1);
        subplot(2,1,1);plot(diff(testDiff))
        subplot(2,1,2);plot(diff(Out.ts))
    end
    
    if negOff < 0 && offset == 0
        Linit.start_time = L.start_time((offset+1):length(L.start_time)) - initT((1):(length(initT)+negOff))'/1000;
        testDiff = L.start_time((offset+1):length(L.start_time)) - initT(1:(length(initT)+negOff))'/1000;
        
        figure(1);
        subplot(2,1,1);plot(diff(testDiff))
        subplot(2,1,2);plot(diff(Out.ts))
        
        
        indexOK = Out.indexOK(1:end-1);
        
        Linit.start_time = Linit.start_time(indexOK);
        D = Out.D;
        D = D(indexOK,:);
    end;
    
    
    %% Get data from Neuralynx
    
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
        ind = find( TTNums == tn )
        fname = [ DataLocation filesep 'TT' num2str(TTNums(ind)) '.ntt']
        tetrodeSpikeTimes = read_cheetah_data(fname);
        load([autoClusterSpikeFiles(TTOrder(ind)).name])
        eval(['tetrodeSet = TT' num2str(TTNums(ind)) '.labels;'])
        
        
        Sc = spikesAuto(Se,session_num,TTNums(ind))
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
        
        experimentType  = 'Classic';
        
        switch experimentType
            case 'Classic'
                % pure correct
                Vis_corr = find( D(:,1) == 1 & D(:,3)==1  );
%                 VisR_corr = find( D(:,1) == 1 & D(:,5)==80 & D(:,2) == -1 );
%                 Vis_corr  = sort( [VisL_corr; VisR_corr] );
                
                % pure correct
                Aud_corr = find( D(:,1) == 1 & D(:,3)==2 );
%                 AudR_corr = find( D(:,1) == 1 & D(:,5)==80 & D(:,2) == -2 );
%                 Aud_corr  = sort( [AudL_corr; AudR_corr] );
                
%                 % mix correct
%                 VisMixL_corr = find( D(:,1) == 1 & D(:,5)==60 & D(:,2) == 1 );
%                 VisMixR_corr = find( D(:,1) == 1 & D(:,5)==60 & D(:,2) == -1 );
%                 VisMix_corr  = sort( [VisMixL_corr; VisMixR_corr] );
%                 
%                 % mix correct
%                 AudMixL_corr = find( D(:,1) == 1 & D(:,5)==60 & D(:,2) == 2 );
%                 AudMixR_corr = find( D(:,1) == 1 & D(:,5)==60 & D(:,2) == -2 );
%                 AudMix_corr  = sort( [AudMixL_corr; AudMixR_corr] );
%                 
%                 % pure incorrect
                 Vis_incorr = find( D(:,1) == 0 & D(:,3)==1 );
%                 VisR_incorr = find( D(:,1) == 0 & D(:,5)==80 & D(:,2) == -1 );
%                 Vis_incorr  = sort( [VisL_incorr; VisR_incorr] );
                
%                 % pure incorrect
                 Aud_incorr = find( D(:,1) == 0 & D(:,3)==2 );
%                 AudR_incorr = find( D(:,1) == 0 & D(:,5)==80 & D(:,2) == -2 );
%                 Aud_incorr  = sort( [AudL_incorr; AudR_incorr] );
                
            case 'New'
                
                % pure correct
                VisL_corr = find( D(:,1) == 1 & D(:,8)==1 & D(:,2) == 1 );
                VisR_corr = find( D(:,1) == 1 & D(:,8)==1 & D(:,2) == -1 );
                Vis_corr  = sort( [VisL_corr; VisR_corr] );
                
                % pure correct
                AudL_corr = find( D(:,1) == 1 & D(:,8)==1 & D(:,2) == 2 );
                AudR_corr = find( D(:,1) == 1 & D(:,8)==1 & D(:,2) == -2 );
                Aud_corr  = sort( [AudL_corr; AudR_corr] );
                
                % mix correct
                VisMixL_corr = find( D(:,1) == 1 & D(:,8)==2 & D(:,2) == 1 );
                VisMixR_corr = find( D(:,1) == 1 & D(:,8)==2 & D(:,2) == -1 );
                VisMix_corr  = sort( [VisMixL_corr; VisMixR_corr] );
                
                % mix correct
                AudMixL_corr = find( D(:,1) == 1 & D(:,8)==2 & D(:,2) == 2 );
                AudMixR_corr = find( D(:,1) == 1 & D(:,8)==2 & D(:,2) == -2 );
                AudMix_corr  = sort( [AudMixL_corr; AudMixR_corr] );
                
                % pure incorrect
                VisL_incorr = find( D(:,1) == 1 & D(:,5)==3 & D(:,2) == 1 );
                VisR_incorr = find( D(:,1) == 1 & D(:,5)==3 & D(:,2) == -1 );
                Vis_incorr  = sort( [VisL_incorr; VisR_incorr] );
                
                % pure incorrect
                AudL_incorr = find( D(:,1) == 1 & D(:,5)==3 & D(:,2) == 2 );
                AudR_incorr = find( D(:,1) == 1 & D(:,5)==3 & D(:,2) == -2 );
                Aud_incorr  = sort( [AudL_incorr; AudR_incorr] );
                
        end;
                
        
        
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
%                 elseif m == 3
%                     f = VisL_corr;
%                     stringVal = 'VisL_corr';
%                     stringVals{m} = stringVal;
%                 elseif m == 4
%                     f = VisR_corr;
%                     stringVal = 'VisR_corr';
%                     stringVals{m} = stringVal;
%                 elseif m == 5
%                     f = AudL_corr;
%                     stringVal = 'AudL_corr';
%                     stringVals{m} = stringVal;
%                 elseif m == 6
%                     f = AudR_corr;
%                     stringVal = 'AudR_corr';
%                     stringVals{m} = stringVal;
%                 elseif m == 7
%                     f = AudMix_corr;
%                     stringVal = 'AudMix_corr';
%                     stringVals{m} = stringVal;
%                 elseif m == 8
%                     f = VisMix_corr;
%                     stringVal = 'VisMix_corr';
%                     stringVals{m} = stringVal;
%                     
                elseif m == 3
                    f = Vis_incorr;
                    stringVal = 'Vis_incorr';
                    stringVals{m} = stringVal;
                elseif m == 4
                    f = Aud_incorr;
                    stringVal = 'Aud_incorr';
                    stringVals{m} = stringVal;
%                 elseif m == 11
%                     f = VisL_incorr;
%                     stringVal = 'VisL_incorr';
%                     stringVals{m} = stringVal;
%                 elseif m == 12
%                     f = VisR_incorr;
%                     stringVal = 'VisR_incorr';
%                     stringVals{m} = stringVal;
%                 elseif m == 13
%                     f = AudL_incorr;
%                     stringVal = 'AudL_incorr';
%                     stringVals{m} = stringVal;
%                 elseif m == 14
%                     f = AudR_incorr;
%                     stringVal = 'AudR_incorr';
%                     stringVals{m} = stringVal;
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
            
            clear Sc Se clcurr_aligncl clcurr_align clcurr_align1 clcurr_align2 clcurr_align3 clcurr_align4 h experimentType

            
        end;
        cd(figpath);
        
        if extt == 0
            save(strcat(Fpath(end-18:end),'_PFCMDSection',num2str(1),'TT',num2str(tn),'.mat'),'-v7.3');
        elseif extt == 1
            save(strcat('PFCMDSection',num2str(1),'TT',num2str(tn),'.mat'),'-v7.3');
        end
        % save(strcat(Fpath(end-18:end),'_PFCMDSection',num2str(cn),'P.mat'),'AudCorr','audCorr','sn');
        clearvars -except sn cn tn numTT MouseID DataLocation behaviorfilename  TT_to_process h experimentType
        close all;
        
    end;
end;

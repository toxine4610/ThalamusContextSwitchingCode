function BuildData_Navi

%%

addpath('C:\Users\Halassalab-CG\Desktop\Session Import');
addpath('C:\Users\Halassalab-CG\Dropbox (Personal)\Rajeev\Rajeev_Code\')
addpath('C:\Users\Halassalab-CG\Dropbox (Personal)\Rajeev\Rajeev_Code\Clustering and Basic Analysis');
addpath('C:\Users\Halassalab-CG\Dropbox (Personal)\Rajeev\Rajeev_Code\Clustering and Basic Analysis\mclust-3.4')
addpath('C:\Users\Halassalab-CG\Dropbox (Personal)\Rajeev\Rajeev_Code\Clustering and Basic Analysis\readcheetahdata')
addpath('C:\Users\Halassalab-CG\Dropbox (Personal)\Rajeev\Rajeev_Code\Clustering and Basic Analysis\readcheetahdata\Optogenetics Toolbox\trunk\@spikesAuto')

%% INIT


clear all;
% MouseID = 'BPM1_3';
DataLocation = 'X:\Ralf\CheetahData\eSSFO1\2018-05-22_16-57-35';
% behaviorfilename.C1 = [ DataLocation filesep 'bpm13_PFCterminal_alltrials_20180518.txt'];
% behaviorfilename.C2 = behaviorfilename.C1;

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


%% CREATE DUMMY SESSION -- NO NEED TO CHANGE
extt = 0; % save to drobo;

session_num =  1;
Se = sessions();
Se = add(Se,'CMO',...
    [2014 04 13 15 11 22],...
    DataLocation,...
    'Neuralynx G',24,...
    'File 6',...
    [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4],... % channels per electrode/home
    [NaN NaN NaN NaN NaN NaN]);

Fpath = Se.folder{session_num};
figpath=strcat(Fpath, filesep, 'GlobalpeakC1');

%% LOAD TTL

%========== Load TTL Values
eventData = read_cheetah_data([Se.folder{session_num} filesep 'Events.nev']);
TTLvalues = unique( eventData.TTLval );

TTLvalues = [0,40];

TTLbehave = 8;     TTLbehave2 = 193;

ct = 0;
for i = 2
    thisTTL = TTLvalues(i);
    L(i-1).start_time = eventData.ts(find((eventData.TTLval == 40 )+(eventData.TTLval == TTLbehave2)))';
end;

%% PULL OUT UNITS
unitCount = 1;

eventData = read_cheetah_data([Se.folder{session_num} '\Events.nev']);

autoClusterSpikeFiles = dir([ DataLocation filesep 'TT*cluster.mat' ] );

for n = 1:length(autoClusterSpikeFiles)
    TTname = autoClusterSpikeFiles(n).name;
    numTT(n) = str2num(TTname(3:(length(TTname)-11)));
end

[TTNums,TTOrder] = sort(numTT);


for i = TTNums
    fname = [ DataLocation filesep 'TT' num2str(TTNums(i)) '.ntt'];
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

%% Manual Cluster

autoclust = 0;


if autoclust == 0
    
    names = dir('*.cut');
    all_tt_nums = [];
    
    
    for n = 1:length(names)
        nameSU = names(n).name;
        nameSU(strfind(nameSU,'.cut'):length(nameSU)) = [];
        nameSU(1:2) = [];
        all_tt_nums = [all_tt_nums str2num(nameSU)];
    end
    
    all_tt_nums = sort(all_tt_nums); %Automatic selection of electrode sets to plot
    
    
    %File Name Adapter
    
    names = dir('*.ntt')
    newNames = names;
    
    for n = 1:length(names)
        newNames(n).name(1) = 'T';
        newNames(n).name(2) = 'T';
        dos(['rename "' names(n).name '" "' newNames(n).name '"']);
    end
    
    names = dir('*.cut')
    newNames = names;
    
    for n = 1:length(names)
        newNames(n).name(1) = 'T';
        newNames(n).name(2) = 'T';
        dos(['rename "' names(n).name '" "' newNames(n).name '"']);
    end
    
%     all_tt_nums = [];
%     
    % == manual import
    fprintf('Extracting Data....');
    
    i = 0;
    
    for n = 1:(length(all_tt_nums))
        curr_tt_num = all_tt_nums(n);
        Sc = spikes(Se,session_num,curr_tt_num);
        is_cluster = 1;
        m = 0;
        while is_cluster == 1
            m = m+1;
            cl_holder = cluster(Sc,m);
            is_full = max(size(cl_holder.timestamp));
            if is_full > 1
                i = i+1;
                eval(sprintf('cl%d = cluster(Sc,m)', i));
                num_seq(i,1:2) = [curr_tt_num m];
                
                
                
            else
                m = m-1;
                is_cluster = 0;
            end
        end
        Sc_unit_count(n,1) = curr_tt_num;
        Sc_unit_count(n,2) = m;
    end
    
    clear cl_holder i curr_tt_num is_cluster is_full m n
    fprintf('....Done!\n');
    
end


clc;
disp('FINISHED');
%%
FigSaveDir = [ DataLocation '\Processed\'];
mkdir(FigSaveDir);

pre = 0.5;
post = 1;



filtWidth = 0.005; bin = 0.010;


for n = 1:size(num_seq,1)
    
    PSTH_Raster_plot =  figure(n);
    set(PSTH_Raster_plot, 'position', [ 2715         564        1742         717])
    set(PSTH_Raster_plot,'color','w');
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultLineLineWidth',1)
    hold on;
    
    figure(n); title(num2str(n));
    
    for i = 1:size(L,2)
        
        clear Ltag f1
        Ltag = L(i);
        f1 = 1:length(Ltag.start_time);
        
        
        eval(['cl' num2str(n) '_align = align(cl' num2str(n) ',Ltag,f1,pre,post);'])
        eval(['clcurr_align = cl' num2str(n) '_align;']);
        SpikeTimesCell =  cl2spikeTimes(clcurr_align);
        
        Spikes(n).SpikeTimes = SpikeTimesCell;
        
        
        [spikeRate,error,errorBounds,spikeMat,time,clCell] = spikeRateEstCltype(clcurr_align,bin,filtWidth,pre,post,Ltag, f1);
        
        
        %         figure(n); subplot(4,2, i + (i-1) ); plotRaster(SpikeTimesCell);xlim([-pre, post]);
        %         figure(n); subplot(4,2, 2*i ); plot( time, spikeRate, 'r');xlim([-pre, post]);
        %         axis normal; box off; set(gca,'tickdir','out','fontsize',16);
        
        figure(n); subplot(1,2, 1 ); plotRaster(SpikeTimesCell);xlim([-pre, post]);
        figure(n); subplot(1,2,2 ); plot( time, spikeRate, 'r');xlim([-pre, post]);
        axis normal; box off; set(gca,'tickdir','out','fontsize',16);
        
        drawnow
        
    end;
    
    print(PSTH_Raster_plot, [FigSaveDir 'PFC_Cell_' num2str(n)], '-dpng');
    close all
    
end;


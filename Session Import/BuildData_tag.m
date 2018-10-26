%%==== process data = optotagging
clear all;
MouseID = 'Constantine';
DataLocation = 'Y:\CheetahData\Constantine\2017-11-15_09-53-47';
% behaviorfilename.C1 = [ DataLocation filesep 'const_20171111_c1_3p.txt'];
% behaviorfilename.C2 = [ DataLocation filesep 'const_20171111_c2_4p_RS.txt'];

cd(DataLocation);
addpath(DataLocation);

switch MouseID
    case 'SOMCRE1'
        TT_to_process = [1:16];
    case 'WT24'
        TT_to_process  = [16];
    case 'Constantine'
        TT_to_process  = [1:18];
end;




%% Make a Temp Session
extt = 0; % save to drobo;

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
    case 'Constantine'
        Se = sessions();
        Se = add(Se,'CMO',...
            [2014 04 13 15 11 22],...
            DataLocation,...
            'autoCluster',32,...
            'File 6',...
            [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4],... % channels per electrode/home
            [NaN NaN NaN NaN NaN NaN]);
end

Fpath = Se.folder{session_num};
figpath=strcat(Fpath,'\GlobalpeakC1');

%%
unitCount = 1;

eventData = read_cheetah_data([Se.folder{session_num} '\Events.nev']);

autoClusterSpikeFiles = dir([ DataLocation filesep 'TT*cluster.mat' ] );

for n = 1:length(autoClusterSpikeFiles)
    TTname = autoClusterSpikeFiles(n).name;
    numTT(n) = str2num(TTname(3:(length(TTname)-11)));
end

[TTNums,TTOrder] = sort(numTT);

for i = [14,15,18]
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

FigSaveDir = [ DataLocation '\Tagged\'];
mkdir(FigSaveDir);

pre = 0.2;
post = 0.4;

% Ltag.start_time = eventData.ts( 955:end ); % 1110
Ltag.start_time = eventData.ts( 840:end-20 ); % 1109
f1 = 1:2:length(Ltag.start_time);
filtWidth = 0.005;

bin = 0.0010;

for n = 1:size(num_seq,1);
    
    eval(['cl' num2str(n) '_align = align(cl' num2str(n) ',Ltag,f1,pre,post);'])
    eval(['clcurr_align = cl' num2str(n) '_align;']);
    SpikeTimesCell =  cl2spikeTimes(clcurr_align);
    
    
    PSTH_Raster_plot =  figure(n);
    set(PSTH_Raster_plot, 'position', [ 2715         564        1742         717])
    set(PSTH_Raster_plot,'color','w');
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultLineLineWidth',1)
    hold on;
    
    [spikeRate,error,errorBounds,spikeMat,time,clCell] = spikeRateEstCltype(clcurr_align,bin,filtWidth,pre,post,Ltag, f1);
    
   figure(n); title(num2str(n));
    figure(n); subplot(1,2,1); plotRaster(SpikeTimesCell);xlim([-0.2, 0.4]);
    figure(n);subplot(1,2,2); plot( time, spikeRate, 'r');xlim([-0.2, 0.4]);axis square; box off; set(gca,'tickdir','out','fontsize',16);
    
    print(PSTH_Raster_plot, [FigSaveDir 'PFC_Cell_' num2str(n+68)], '-dpng');
    close all
    
end;





function [ initT, trial, laser, LEDTBAud, right_trials, left_trials, stim_lengths, catchh, audioStim, visStim, audioStimDistract, visStimDistract, corr, inc, ts, cueL, cueR, ITI, LEDstate] = parseText( filename, MouseID)

trial  = [];
laser  = [];
double_laser = [];
right_trials  = [];
left_trials   = [];
catchh = [];
corr   = [];
inc    = [];
ts     = [];
charge = [];
charge_time = [];
stim_lengths = [];
stim_duration = [];
cueL =[];
cueR = [];
auto_ind = [];
LEDTGVis = [];
LEDTBAud = [];
audioStim = [];
visStim = [];
audioStimDistract = [];
visStimDistract = [];

state     = 'sleep';
lastTs    = 0;
currentTs = 0;
trialInd  = 0;
stimtime  = 0;

fname = filename;
[behavior] = textread(fname,'%s',-1,'delimiter','\t');



%%
switch MouseID
    case 'SOMCRE1'
       LEDstate = behavior(2:2:end);
        timestamps = behavior(1:2:end-1);
    case 'Constantine'
         LEDstate = behavior(2:2:end);
        timestamps = behavior(1:2:end-1);
    case 'WT24'
        LEDstate = behavior(2:2:end);
        timestamps = behavior(1:2:end-1);
    case 'BPM1_3'
        timestamps = behavior(1:2:end-1);
        LEDstate = behavior(2:2:end);
end
% [led,ts] = cleanBehavior(behavior);



% %separating file into LED state and timestamps
%
%  LEDstate = led';
%
% %convert to an aray of doubles
%
%  timestamps   = ts';



%separating file into LED state and timestamps
% LEDstate = behavior(1:2:end-1)
%
% %convert to an aray of doubles
% timestamps = behavior(2:2:end)

%%

Timestamp = [];
length(Timestamp)

gapfind = 0;
gapOff = 0;

for i = 1:length(timestamps);
    timestamp = cell2mat(timestamps(i,1));
    timestamp = int32(str2num(timestamp));
    timestamp = double(timestamp);
    if (i > 2 && gapfind < 1)
        if timestamp < Timestamp(i-1)
            gapOff = Timestamp(i-1);
            gapfind = gapfind+1;
        end
    end
    timestamp = timestamp+gapOff;
    Timestamp = [Timestamp;timestamp];
end


for ii = 1:length(LEDstate)
    transition = LEDstate{ii};
    
    currentTs  = Timestamp(ii);
    
    switch state
        case 'init'
            if strfind(transition, 'yellow laser on')
                laser = [laser trialInd];
            end
            if strcmp(transition, 'both laser on')
                double_laser = [double_laser trialInd];
            end
            if strfind(transition, 'autostart')
                state = 'trial';
            end
        case 'trial'
            stimtime = currentTs;
            if strfind(transition,'cueL')
                cueL = [cueL trialInd];
                state = 'trial';
            elseif strfind(transition, 'cueR')
                cueR = [cueR trialInd];
                state = 'trial';
            else
                state = 'trial';
            end
            if strfind(transition, 'stim l on') == 1
                stim_duration = transition(11:13);
                stim_duration = int32(str2num(stim_duration));
                stim_duration = double(stim_duration);
                left_trials = [left_trials trialInd];
                ts   = [ts currentTs];
                state = 'stim type';
            elseif strfind(transition, 'stim r on') == 1
                stim_duration = transition(11:13);
                stim_duration = int32(str2num(stim_duration));
                stim_duration = double(stim_duration);
                right_trials = [right_trials trialInd];
                ts   = [ts currentTs];
                state = 'stim type';
            elseif strfind(transition, 'stim left on')
                left_trials = [left_trials trialInd];
                ts   = [ts currentTs];
                state = 'stim type';
            elseif strfind(transition, 'stim right on')
                right_trials = [right_trials trialInd];
                ts   = [ts currentTs];
                state = 'stim type';
            elseif strfind(transition, 'catch')
                catchh = [catchh, trialInd];
                ts   = [ts currentTs];
                state = 'stim type';
            end
            stim_lengths = [stim_lengths stim_duration];
        case 'stim type'
            if strcmp(transition,'audio stim')
                audioStim = [audioStim trialInd];
                state = 'choice';
            elseif strcmp(transition,'visual stim')
                visStim = [visStim trialInd];
                state = 'choice';
            elseif strfind(transition,'audio stim')&strfind(transition,'distractor')
                audioStimDistract = [audioStimDistract trialInd];
                state = 'choice';
            elseif strfind(transition,'visual stim')&strfind(transition,'distractor')
                visStimDistract = [visStimDistract trialInd];
                state = 'choice';
            end
            latency = currentTs - stimtime;
            if strfind(transition, 'reward delivered');
                trial = [trial latency];
                corr = [corr trialInd];
                state = 'sleep';
            elseif strfind(transition, 'stim off')
                trial = [trial latency];
                inc = [inc trialInd];
                state = 'sleep';
            end
        case 'choice'
            latency = currentTs - stimtime;
            if strfind(transition, 'reward delivered');
                trial = [trial latency];
                corr = [corr trialInd];
                state = 'sleep';
            elseif strfind(transition, 'stim off')
                trial = [trial latency];
                inc = [inc trialInd];
                state = 'sleep';
            end
        case 'sleep'
            if strfind(transition, 'LEDTGWhite On')
                LEDTGVis = [LEDTGVis trialInd];
                state = 'init';
                trialInd = trialInd + 1;
            elseif strfind(transition, 'LEDTG On')
                LEDTGVis = [LEDTGVis trialInd];
                state = 'init';
                trialInd = trialInd + 1;
            elseif strfind(transition, 'LEDTB On')
                LEDTBAud = [LEDTBAud trialInd];
                state = 'init';
                trialInd = trialInd + 1;
            elseif strfind(transition, 'initiation Visual')
                LEDTGVis = [LEDTGVis trialInd];
                state = 'init';
                trialInd = trialInd + 1;
            elseif strfind(transition, 'initiation Audio')
                LEDTBAud = [LEDTBAud trialInd];
                state = 'init';
                trialInd = trialInd + 1;
            end
            
    end;
    
    %     if strfind(transition,'repeat cue different')
    % %         disp("here");
    %         cue_trials_different(trialInd) = 1;
    %     end
    %
    %      if strfind(transition,'repeat cue same')
    % %         disp("here");
    %         cue_trials_same(trialInd) = 1;
    %     end
    
    lastTs = currentTs;
    
end


%%
auto_ind = Timestamp(find(strcmp(LEDstate, 'autostart')));
initInd = Timestamp(find(strcmp(LEDstate, 'initiation')));

if length(initInd) < 1
        initIndTB = Timestamp(find(strcmp(LEDstate, 'LEDTB On')));
        initIndTG = Timestamp(find(strcmp(LEDstate, 'LEDTG On')));
    initIndTB = Timestamp(find(~cellfun(@isempty,(strfind(LEDstate, 'LEDTB On')))));
    initIndTG = Timestamp(find(~cellfun(@isempty,(strfind(LEDstate, 'LEDTG On')))));
    initInd = sortrows([initIndTB;initIndTG]);
end
length(initIndTB)
for n = 1:length(auto_ind)
    n
    initIndHold = auto_ind(n) - initInd
    initIndHold = initIndHold(find(initIndHold > 0))
    initT(n) = min(initIndHold);
end

%  initT = computeTimeStamps(LEDstate, timestamps);
 
 m = mode(round(initT));
 ind_common = find( initT <= m+60 & initT >= m-60 );
 ind_uncommon = setdiff( 1:length(initT), ind_common );
 if ~isempty(ind_uncommon)
     for i = 1:length(ind_uncommon)
         this_ = initT(ind_uncommon(i));
         diff  = m - this_;
         initT( ind_uncommon(i) ) = m;
     end;
 end;
 

length(initT)
auto_ind = find(strcmp(LEDstate, 'autostart'));
for n = 2:length(auto_ind);
    ITIval = Timestamp(auto_ind(n)) - Timestamp(auto_ind(n - 1));
    ITI(n) = ITIval;
end

length(laser)
laser = unique(laser);
double_laser = unique(double_laser);





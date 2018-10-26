function [initT, trial, right_trials, left_trials,...
    audioStim, audioStimDistract, visStim, visStimDistract,...
    corr, inc, missed, ts, ITI, LEDstate, ErrInLeft, ErrInRight, ErrOutLeft, ErrOutRight, holdthroughtrials, accumulationtrials, mixedTrials, mixedPercent] =...
    ParseEngineCrossModal4Poke_TreeShrewNoLaser2(filename);

%%
trial  = [];
right_trials  = [];
left_trials   = [];
corr   = [];
inc    = [];
ts     = [];
audioStim = [];
audioStimDistract = [];
visStim = [];
visStimDistract = [];
missed = [];
ErrInLeft = [];
ErrInRight = []; 
ErrOutLeft = []; 
ErrOutRight = [];
trialtype = []; %0=click, 1=holdthrough
mixedTrials = [];
mixedPercent = [];


state     = 'sleep';
lastTs    = 0;
currentTs = 0;
trialInd  = 0;
stimtime  = 0;

fname = filename;
[behavior] = textread(fname,'%s',-1,'delimiter','\t');

%separating file into LED state and timestamps
LEDstate = behavior(1:2:end-1);

%convert to an aray of doubles
timestamps = behavior(2:2:end);

%%
Timestamp = [];

gapfind = 0;
gapOff = 0;

for i = 1:length(timestamps)
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
    
    if length(transition) == 16 && strfind(transition, 'mixed trial') && ~ismember(trialInd,mixedTrials)
        mixedTrials = [mixedTrials trialInd];
        mixedPercent = [mixedPercent str2double(transition(13:end))];
    end
            
    switch state
        case 'init'
            if strfind(transition,'LEDTB On')
                audioStim = [audioStim trialInd];
            elseif strfind(transition, 'LEDTG On')
                visStim = [visStim trialInd];
            elseif strfind(transition, 'Hold Through Trial')
                trialtype(trialInd) = 1;
            elseif strfind(transition, 'Click Trial')
                trialtype(trialInd) = 0;
            elseif strfind(transition, 'autostart')
                state = 'trial';
            end
        case 'trial'
            stimtime = currentTs;           
            if strfind(transition, 'stim left on') == 1
               left_trials = [left_trials trialInd];
               ts   = [ts currentTs];
               state = 'stim type';
            elseif strfind(transition, 'stim right on') == 1
               right_trials = [right_trials trialInd];
               ts   = [ts currentTs];
               state = 'stim type';
            end
        case 'stim type'
        
            if strcmp(transition,'audio stim')
                audioStim = [audioStim trialInd];
                state = 'choice';
            elseif strcmp(transition,'audio stim distractor')
                audioStimDistract = [audioStimDistract trialInd];
                state = 'choice';
            elseif strcmp(transition,'visual stim')
                visStim = [visStim trialInd];
                state = 'choice';
            elseif strcmp(transition,'visual stim distractor')
                visStimDistract = [visStimDistract trialInd];
                state = 'choice';
            end
            latency = currentTs - stimtime;  
            if strfind(transition, 'reward delivered')
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
            if strfind(transition, 'reward delivered')
                trial = [trial latency];
                corr = [corr trialInd];
                state = 'sleep';
            elseif strfind(transition, 'stim off')
                trial = [trial latency];
                inc = [inc trialInd];
                state = 'sleep';
            elseif strfind(transition, 'wrong poke!')
                trial = [trial latency];
                inc = [inc trialInd];
                if strfind(transition, 'Outer left poked')
                    ErrOutLeft = [ErrOutLeft trialInd];
                elseif strfind(transition, 'Outer right poked')
                    ErrOutRight = [ErrOutRight trialInd];
                elseif strfind(transition, 'Center left poked')
                    ErrInLeft = [ErrInLeft trialInd];
                elseif strfind(transition, 'Center right poked')
                    ErrInRight = [ErrInRight trialInd];
                end
                state = 'sleep';
            elseif strfind(transition, 'missed your shot!')
                trial = [trial latency];
                missed = [missed trialInd];
                state = 'sleep';
            end
        case 'sleep'  
            if strfind(transition, 'trial available')
                trialInd = trialInd + 1;
            elseif strfind(transition, 'LEDTGWhite On')
                state = 'init';
            end  
    end
    lastTs = currentTs;
end

auto_ind = find(strcmp(LEDstate, 'autostart') );
initT = Timestamp(auto_ind) - Timestamp(auto_ind - 1);

for n = 2:length(auto_ind)
    ITIval = Timestamp(auto_ind(n)) - Timestamp(auto_ind(n - 1));
    ITI(n) = ITIval; 
end

holdthroughtrials = find(trialtype == 1);
accumulationtrials = find(trialtype == 0);
audioStim = unique(audioStim);
visStim = unique(visStim);

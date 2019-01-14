
function [glmtrial, unitOfTime, binSize, nTrials, binfun, numMD, numPFC, RejMD, RejPFC, CMIPFC, CMIMD] = GLMpreprocess(dataName, behavName)

dataRep = '/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/RuleSwitchDataBase/'

load( [dataRep dataName] );
load( [dataRep behavName] );


addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/RuleSwitchDataBase/WT19');
saveFlag = 0;

%% Combine two contexts together
% Z(:,9) = Context (1 if 1 and 2 if 2)

Z_Context1  = [ Z_R1C1 ; Z_R2C1 ];
Z_Context2  = [ Z_R1C2 ; Z_R2C2 ];
Z_Context1 = [Z_Context1 ones(size(Z_Context1,1),1)];
Z_Context2 = [Z_Context2 2*ones(size(Z_Context2,1),1)];

% Z contains both contexts
Z = [Z_Context1; Z_Context2];

%% PFC
numPFC = length(Spfc_singlesession);
pfc_spiketimes = struct();

ct = 0; RejPFC = [];

for i = 1:numPFC
    
    SpikeTimes = cell(size(Z,1),1);
    
    correct = 1; rule = 1; context = 1;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Spfc_singlesession(i).SpikeTimes_R1C1;
    
    correct = 1; rule = 2; context = 1;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Spfc_singlesession(i).SpikeTimes_R2C1;
    
    correct = 1; rule = 1; context = 2;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Spfc_singlesession(i).SpikeTimes_R1C2;
    
    correct = 1; rule = 2; context = 2;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Spfc_singlesession(i).SpikeTimes_R2C2;
    
    f4 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Spfc_singlesession(i).SpikeTimes_R2C2 ) ~= 0 ) )./numel(Spfc_singlesession(i).SpikeTimes_R2C2);
    f3 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Spfc_singlesession(i).SpikeTimes_R1C2 ) ~= 0 ) )./numel(Spfc_singlesession(i).SpikeTimes_R1C2);
    
    f2 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Spfc_singlesession(i).SpikeTimes_R2C1 ) ~= 0 ) )./numel(Spfc_singlesession(i).SpikeTimes_R2C1);
    f1 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Spfc_singlesession(i).SpikeTimes_R1C1 ) ~= 0 ) )./numel(Spfc_singlesession(i).SpikeTimes_R1C1);
    
    if length( find( [f1,f2,f3,f4] <= 0.25 ) ) == 0
        ct = ct+1;
        
        
        
        
        %     correct = 0; rule = 1; context = 1;
        %     trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
        %     SpikeTimes(trialSubset) = Spfc_singlesession_inc(i).SpikeTimes_R1C1;
        %
        %     correct = 0; rule = 2; context = 1;
        %     trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
        %     SpikeTimes(trialSubset) = Spfc_singlesession_inc(i).SpikeTimes_R2C1;
        %
        %     correct = 0; rule = 1; context = 2;
        %     trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
        %     SpikeTimes(trialSubset) = Spfc_singlesession_inc(i).SpikeTimes_R1C2;
        %
        %     correct = 0; rule = 2; context = 2;
        %     trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
        %     SpikeTimes(trialSubset) = Spfc_singlesession_inc(i).SpikeTimes_R2C2;
        %
        pfc_spiketimes(ct).SpikeTimes = SpikeTimes;
        pfc_spiketimes(ct).filename_correct = Spfc_singlesession(i).filename;
        %     pfc_spiketimes(i).filename_incorrect = Spfc_singlesession_inc(i).filename;
        pfc_spiketimes(ct).unitNbr_correct = Spfc_singlesession(i).unitNbr;
        %     pfc_spiketimes(i).unitNbr_incorrect = Spfc_singlesession_inc(i).unitNbr;
        %     pfc_spiketimes(i).tetrodeNumber_incorrect = Spfc_singlesession_inc(i).ttNbr;
    else
        RejPFC = [RejPFC,i];
    end;
end

numPFC = ct;

%%
numMD = length(Smd_singlesession);
md_spiketimes = struct();

ct = 0; RejMD = [];

for i = 1:numMD
    
    SpikeTimes = cell(size(Z,1),1);
    
    correct = 1; rule = 1; context = 1;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Smd_singlesession(i).SpikeTimes_R1C1;
    
    correct = 1; rule = 2; context = 1;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Smd_singlesession(i).SpikeTimes_R2C1;
    
    correct = 1; rule = 1; context = 2;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Smd_singlesession(i).SpikeTimes_R1C2;
    
    correct = 1; rule = 2; context = 2;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Smd_singlesession(i).SpikeTimes_R2C2;
    
    f4 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Smd_singlesession(i).SpikeTimes_R2C2 ) ~= 0 ) )./numel(Smd_singlesession(i).SpikeTimes_R2C2);
    f3 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Smd_singlesession(i).SpikeTimes_R1C2 ) ~= 0 ) )./numel(Smd_singlesession(i).SpikeTimes_R1C2);
    
    f2 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Smd_singlesession(i).SpikeTimes_R2C1 ) ~= 0 ) )./numel(Smd_singlesession(i).SpikeTimes_R2C1);
    f1 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Smd_singlesession(i).SpikeTimes_R1C1 ) ~= 0 ) )./numel(Smd_singlesession(i).SpikeTimes_R1C1);
    
    if length( find( [f1,f2,f3,f4] <= 0.25 ) ) == 0
        ct = ct+1;
        
        %     correct = 0; rule = 1; context = 1;
        %     trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
        %     SpikeTimes(trialSubset) = Smd_singlesession_inc(i).SpikeTimes_R1C1;
        %
        %     correct = 0; rule = 2; context = 1;
        %     trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
        %     SpikeTimes(trialSubset) = Smd_singlesession_inc(i).SpikeTimes_R2C1;
        %
        %     correct = 0; rule = 1; context = 2;
        %     trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
        %     SpikeTimes(trialSubset) = Smd_singlesession_inc(i).SpikeTimes_R1C2;
        %
        %     correct = 0; rule = 2; context = 2;
        %     trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
        %     SpikeTimes(trialSubset) = Smd_singlesession_inc(i).SpikeTimes_R2C2;
        %
        md_spiketimes(ct).SpikeTimes = SpikeTimes;
        md_spiketimes(ct).filename_correct = Smd_singlesession(i).filename;
        %     md_spiketimes(i).filename_incorrect = Smd_singlesession_inc(i).filename;
        md_spiketimes(ct).unitNbr_correct = Smd_singlesession(i).unitNbr;
        %     md_spiketimes(i).unitNbr_incorrect = Smd_singlesession_inc(i).unitNbr;
        %     md_spiketimes(i).tetrodeNumber_incorrect = Smd_singlesession_inc(i).ttNbr;
    else
        RejMD = [RejMD,i];
    end;
end

numMD = ct;

%% Combine spikes and behavior into a trial structure - all in MS
binSize = 1; %ms
glmtrial = struct();
unitOfTime = 'ms';
nTrials = size(Z,1);
binfun = @(t)(t==0)+ceil(t/binSize);

% loop over trials, bin behavioral things into correct bins and add spike
% trains correctly
for i = 1:nTrials
    
    % trial start -> Z(i,4)
    trial_start = Z(i,4);
    
    % trial duration = Z(i,7);
    
    % NEED TO CHECK WITH RAJEEV ABOUT THIS
    %     % trial duration in seconds = Z(i,4)/10000 + Z(i,7)/1000
    %     duration = Z(i,4)/10000 + Z(i,7)/1000;
    
    % trial duration in ms= Z(i,4)/10 + Z(i,7)
    %     trial_start_ms = -Z(i,4)/10;
    trial_start_ms = 0; % start analysis at -5-0 ms
    %     trial_end_ms = Z(i,3);
    trial_end_ms = Z(i,7);
    %     trial_end_ms = 700;
    duration_ms = trial_end_ms - trial_start_ms;
    glmtrial(i).duration = duration_ms;
    
    % context
    glmtrial(i).context = Z(i,9);
    
    glmtrial(i).R1C1 = [];
    glmtrial(i).R2C1 = [];
    glmtrial(i).R1C2 = [];
    glmtrial(i).R2C2 = [];
    
    % cue onset and offset - cue onset always at 0 ms and cue offset always
    % at 100 ms
    cueon = (0-trial_start_ms); % ms
    cueoff = (100 - trial_start_ms); % ms
    
    glmtrial(i).cueon = cueon; % binfun(cueon);
    glmtrial(i).cueoff = cueoff; % binfun(cueoff);
    
    % rule - Z(:, 2) = 1 - vision, 2 - audition
    glmtrial(i).vision = [];
    glmtrial(i).audition = [];
    if Z(i,2)==1 % attend to vision rule
        glmtrial(i).vision = 0;
    elseif Z(i,2)==2 % attend to audition rule
        glmtrial(i).audition = 0;
    end
    %     glmtrial(i).rule = Z(i,2); % 1 = vision, 2 = audition
    
    % cue - low-pass or high-pass noise cue
    % do this like rule -> need to infer cue from rule and context
    % in "rule-switch" condition [R,C]
    % -> context 1, vision = high-pass noise [1,1]
    % -> context 1, audition = low-pass noise [2,1]
    % -> context 2, vision = vision = low-pass noise [1,2]
    % -> context 2, audition = high-pass noise [2,2]
    glmtrial(i).lowpasscue = [];
    glmtrial(i).highpasscue = [];
    if Z(i,9)==1 % context 1
        if Z(i,2)==1 % vision
            %             glmtrial(i).type = 11;
            glmtrial(i).R1C1 = cueon;
            glmtrial(i).highpasscue = 0;
        elseif Z(i,2)==2 % audition
            %             glmtrial(i).type = 21;
            glmtrial(i).R2C1 = cueon;
            glmtrial(i).lowpasscue = 0;
        end
    elseif Z(i,9)==2 % context 2
        if Z(i,2)==1 % vision
            %             glmtrial(i).type = 12;
            glmtrial(i).R1C2 = cueon;
            glmtrial(i).lowpasscue = 0;
        elseif Z(i,2)==2 % audition
            %             glmtrial(i).type = 22;
            glmtrial(i).R2C2 = cueon;
            glmtrial(i).highpasscue = 0;
        end
    end
    
    % reward
    glmtrial(i).reward = Z(i,5);
    
    % choice time = end of trial, which is trial duration
    %     glmtrial(i).choice =  Z(i,7);
    
    % spike trains
    % for each neuron, PFC and MD - get spike trains
    for n = 1:numPFC
        temp_spiketimes = pfc_spiketimes(n).SpikeTimes{i}*1000 - trial_start_ms;
        glmtrial(i).(['PFCUnit' num2str(n)]) = temp_spiketimes(temp_spiketimes>=0 & temp_spiketimes<glmtrial(i).duration);
    end
    
    for n = 1:numMD
        temp_spiketimes = md_spiketimes(n).SpikeTimes{i}*1000 - trial_start_ms;
        glmtrial(i).(['MDUnit' num2str(n)]) = temp_spiketimes(temp_spiketimes>=0 & temp_spiketimes<glmtrial(i).duration);
    end
    
end

%% compute CMI

keptPFC = setdiff( 1:numel(Spfc_singlesession), RejPFC);
keptMD  = setdiff( 1:numel(Smd_singlesession), RejMD);

range = [ 0.0, 0.6 ];
bin = 0.0010;
filtWidth = 0.08;
CMI = @(X,Y) (X-Y)./(X+Y);

for i = 1:numel(Spfc_singlesession)
    clear C1 C2
    C1 = [Spfc_singlesession(i).SpikeTimes_R1C1, Spfc_singlesession(i).SpikeTimes_R2C1];
    
    [C1, ~, time, ~] = makeSpikeRates(C1, range , bin, filtWidth);
    
    C2 = [Spfc_singlesession(i).SpikeTimes_R1C2, Spfc_singlesession(i).SpikeTimes_R2C2];
    
    [C2, ~, time, ~] = makeSpikeRates(C2, range , bin, filtWidth);
    
    CMIPFC(i) = CMI( nanmean(C1), nanmean(C2) );
end;

for i = 1:numel(Smd_singlesession)
    clear C1 C2
    C1 = [Smd_singlesession(i).SpikeTimes_R1C1, Smd_singlesession(i).SpikeTimes_R2C1];
    
    [C1, ~, time, ~] = makeSpikeRates(C1, range , bin, filtWidth);
    
    C2 = [Smd_singlesession(i).SpikeTimes_R1C2, Smd_singlesession(i).SpikeTimes_R2C2];
    
    [C2, ~, time, ~] = makeSpikeRates(C2, range , bin, filtWidth);
    
    CMIMD(i) = CMI( nanmean(C1), nanmean(C2) );
end;

CMIPFC = CMIPFC( keptPFC );
CMIMD  = CMIMD( keptMD );




%%
if saveFlag == 1
    save('RajeevDataGLM_Session2','glmtrial','unitOfTime','binSize','nTrials','binfun')
end;

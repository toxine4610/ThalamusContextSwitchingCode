% single session - WT19

addpath('/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/RuleSwitchDataBase/WT19');

%% Load behavior
% Z_Context1, Z_Context2, Z_R1C1, Z_R1C2, Z_R2C1, Z_R2C2
% Z contains trial by trial data for context 1, context 2, or correct
% trials in context 1 (R1C1, R2C2) or in context 2 (R1C2, R2C2)
% Z(:, 1) = 0 - incorrect; 1- correct
% Z(:, 2) = 1 - vision, 2 - audition
% Z(:, 3) = hold duration
% Z(:, 4) =  time to initiate after broadband noise (t_start of broad band
%           noise - t_initiation)
% Z(:, 5)  = 1 - collect reward on R, 0 - collect reward on L
% Z(:, 6) = IGNORE 
% Z(:, 7) =  total trial time in ms from initiation
% Z(:, 8) = smoothed performance
path = '/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/RuleSwitchDataBase/WT19/GLMSharedFolder/'
load([ path 'WT19_SingleSession_Behavior.mat'])

%% Combine two contexts together 
% Z(:,9) = Context (1 if 1 and 2 if 2)
Z_Context1  = [ Z_R1C1 ; Z_R2C1 ];
Z_Context2  = [ Z_R1C2 ; Z_R2C2 ];
Z_Context1 = [Z_Context1 ones(size(Z_Context1,1),1)];
Z_Context2 = [Z_Context2 2*ones(size(Z_Context2,1),1)];

% Z contains both contexts
Z = [Z_Context1; Z_Context2];

%% Load spikes
% Smd_singlesession, Spfc_singlesession, Smd_singlesession_inc, Spfc_singlesession_inc
% Spike data format
% S[brain area]_singleSession[_inc] -> are structures - (i) is the i^th
% neuron. Have 25 MD and 40 PFC
% e.g. Spfc_singlesession(i) is the responses of the i^th neuron during
% correct trials. The fields are the filename, unit number, tetrode (?)
% number, and spike times on 4 different classes of trials - a combination
% of rule 1 or 2 and context 1 or 2.
load([path 'WT19_SingleSession_RuleSwitched_SpikeTimes.mat'])
load([path 'WT19_SingleSession_Incorrect_SpikeTimes.mat'])


%% Combine spikes
% want two spike structures, one PFC and one MD
% want trials in order of Z (behavior). 
% rule 1 is vision, rule 2 is audition
% Z(:,2)==1 & Z(:,9) == 1 -> Rule 1, Context 1
% Z(:,2)==1 & Z(:,9) == 2 -> Rule 1, Context 2
% Z(:,2)==2 & Z(:,9) == 1 -> Rule 2, Context 1
% Z(:,2)==2 & Z(:,9) == 2 -> Rule 2, Context 2

% need -> low-pass noise + vision 
% need -> low-pass noise + audition
% need -> high-pass noise + vision
% need -> high-pass noise + audition
% R1C1 ---> Cue = Blue noise; Rule = attend to vision ( in context 1)
% R2C1 ---> Cue = Brown noise; Rule = attend to audition ( in context 1)
% R1C2 ---> Cue = Brown noise; Rule = attend to vision ( in context 2)
% R2C2 ---> Cue = Blue noise; Rule = attend to audition ( in context 2)

% Z(:,1) == 0 means incorrect
% Z(:,1) == 1 means correct

%% PFC
numPFC = length(Spfc_singlesession);
pfc_spiketimes = struct();

ct = 0;

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
    end;
end

%%
numMD = length(Smd_singlesession);
md_spiketimes = struct();

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
    md_spiketimes(i).SpikeTimes = SpikeTimes;
    md_spiketimes(i).filename_correct = Smd_singlesession(i).filename;
%     md_spiketimes(i).filename_incorrect = Smd_singlesession_inc(i).filename;
    md_spiketimes(i).unitNbr_correct = Smd_singlesession(i).unitNbr;
%     md_spiketimes(i).unitNbr_incorrect = Smd_singlesession_inc(i).unitNbr;
%     md_spiketimes(i).tetrodeNumber_incorrect = Smd_singlesession_inc(i).ttNbr;
    
end
            
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
%%
save('RajeevDataGLM_Session2','glmtrial','unitOfTime','binSize','nTrials','binfun')





function [glmtrial, unitOfTime, binSize, nTrials, binfun, numMD, numPFC, RejMD, RejPFC, keptPFC, keptMD, CMIPFC, CMIMD, Z] = packageData_for_glm(dataName)

if ~ismac
    dataRep = 'C:\Users\Halassalab-CG\Dropbox\Rajeev\ForContextSwitchProject\DoubleCueDatabase\BPM1_3\UniLatMDHalo-Low\';
elseif ismac
%     dataRep = '/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/DoubleCueDatabase/SOMCre/DualModalityCue/';
    dataRep = '/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/DoubleCueDatabase/BPM1_3/ContextSwitchHalo/UniLatMD-Low';
    
    

end

load( [dataRep dataName] );

%% ==== Package Trials
Z_C1(:,10) = 0;

indCorrVis = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 1 );
indCorrAud = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 2 );

ZVis = Z_C1(indCorrVis, :);
ZAud = Z_C1(indCorrAud, :);

indR1C1 = find( ZVis(:,9) == 0);
indR2C1 = find( ZAud(:, 9) == 0);

indR1C2 = find( ZVis(:,9) == 3 );
indR2C2 = find( ZAud(:, 9) == 3 );


for i = 1:numel(Spfc)
    Spfc_singlesession(i).unitNbr         = i;
    Spfc_singlesession(i).SpikeTimes_R1C1 = Spfc(i).SpikeTimes_R1C1(indR1C1);
    Spfc_singlesession(i).SpikeTimes_R1C2 = Spfc(i).SpikeTimes_R1C1(indR1C2);
    
    Spfc_singlesession(i).SpikeTimes_R2C1 = Spfc(i).SpikeTimes_R2C1(indR2C1);
    Spfc_singlesession(i).SpikeTimes_R2C2 = Spfc(i).SpikeTimes_R2C1(indR2C2);
    
end;


for i = 1:numel(Smd)
    Smd_singlesession(i).unitNbr         = i;
    Smd_singlesession(i).SpikeTimes_R1C1 = Smd(i).SpikeTimes_R1C1(indR1C1);
    Smd_singlesession(i).SpikeTimes_R1C2 = Smd(i).SpikeTimes_R1C1(indR1C2);
    
    Smd_singlesession(i).SpikeTimes_R2C1 = Smd(i).SpikeTimes_R2C1(indR2C1);
    Smd_singlesession(i).SpikeTimes_R2C2 = Smd(i).SpikeTimes_R2C1(indR2C2);
    
end;

%% ===== Organize Z matrix

indCorr = find(Z_C1(:,1) == 1);
Z       = Z_C1(indCorr,:);
% indLas  = find(Z(:,10) == 1);
% Z       = Z(indLas,:);


%% PFC ====================================================================

numPFC = length(Spfc_singlesession);
pfc_spiketimes = struct();

ct = 0; RejPFC = [];

for i = 1:numPFC
    
    SpikeTimes = cell(size(Z,1),1);
    
    correct = 1; rule = 1; context = 0;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Spfc_singlesession(i).SpikeTimes_R1C1;
    
    correct = 1; rule = 2; context = 0;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Spfc_singlesession(i).SpikeTimes_R2C1;
    
    correct = 1; rule = 1; context = 3;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Spfc_singlesession(i).SpikeTimes_R1C2;
    
    correct = 1; rule = 2; context = 3;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Spfc_singlesession(i).SpikeTimes_R2C2;
    
    f4 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Spfc_singlesession(i).SpikeTimes_R2C2 ) ~= 0 ) )./numel(Spfc_singlesession(i).SpikeTimes_R2C2);
    f3 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Spfc_singlesession(i).SpikeTimes_R1C2 ) ~= 0 ) )./numel(Spfc_singlesession(i).SpikeTimes_R1C2);
    
    f2 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Spfc_singlesession(i).SpikeTimes_R2C1 ) ~= 0 ) )./numel(Spfc_singlesession(i).SpikeTimes_R2C1);
    f1 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Spfc_singlesession(i).SpikeTimes_R1C1 ) ~= 0 ) )./numel(Spfc_singlesession(i).SpikeTimes_R1C1);
    
    if length( find( [f1,f2,f3,f4] <= 0.25 ) ) == 0
        ct = ct+1;
        pfc_spiketimes(ct).SpikeTimes = SpikeTimes;
        %         pfc_spiketimes(ct).filename_correct = Spfc_singlesession(i).filename;
        pfc_spiketimes(ct).unitNbr_correct = Spfc_singlesession(i).unitNbr;
        
    else
        RejPFC = [RejPFC,i];
    end;
end

numPFC = ct;

%% MD =====================================================================
numMD = length(Smd_singlesession);
md_spiketimes = struct();

ct = 0; RejMD = [];

for i = 1:numMD
    
    SpikeTimes = cell(size(Z,1),1);
    
    correct = 1; rule = 1; context = 0;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Smd_singlesession(i).SpikeTimes_R1C1;
    
    correct = 1; rule = 2; context = 0;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Smd_singlesession(i).SpikeTimes_R2C1;
    
    correct = 1; rule = 1; context = 3;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Smd_singlesession(i).SpikeTimes_R1C2;
    
    correct = 1; rule = 2; context = 3;
    trialSubset = Z(:,1)==correct & Z(:,2) == rule & Z(:,9) == context;
    SpikeTimes(trialSubset) = Smd_singlesession(i).SpikeTimes_R2C2;
    
    f4 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Smd_singlesession(i).SpikeTimes_R2C2 ) ~= 0 ) )./numel(Smd_singlesession(i).SpikeTimes_R2C2);
    f3 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Smd_singlesession(i).SpikeTimes_R1C2 ) ~= 0 ) )./numel(Smd_singlesession(i).SpikeTimes_R1C2);
    
    f2 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Smd_singlesession(i).SpikeTimes_R2C1 ) ~= 0 ) )./numel(Smd_singlesession(i).SpikeTimes_R2C1);
    f1 = length( find( cellfun( @(x) length( find( x >=0 & x <= 0.5 ) ), Smd_singlesession(i).SpikeTimes_R1C1 ) ~= 0 ) )./numel(Smd_singlesession(i).SpikeTimes_R1C1);
    
    if length( find( [f1,f2,f3,f4] <= 0.25 ) ) == 0
        ct = ct+1;
        
        md_spiketimes(ct).SpikeTimes = SpikeTimes;
        %         md_spiketimes(ct).filename_correct = Smd_singlesession(i).filename;
        md_spiketimes(ct).unitNbr_correct = Smd_singlesession(i).unitNbr;
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

for i = 1:nTrials
    
    trial_start_ms = 0;
    trial_end_ms = 800;
    duration_ms = trial_end_ms - trial_start_ms;
    glmtrial(i).duration = duration_ms;
    
    glmtrial(i).context = Z(i,9);
    
    glmtrial(i).R1C1 = [];
    glmtrial(i).R2C1 = [];
    glmtrial(i).R1C2 = [];
    glmtrial(i).R2C2 = [];
    
    % cue onset and offset - cue onset always at 0 ms and cue offset always
    % at 100 ms
    cueon = (50); % ms
    cueoff = (250); % ms
    
    glmtrial(i).cueon = cueon; % binfun(cueon);
    glmtrial(i).cueoff = cueoff; % binfun(cueoff);
    
    
    glmtrial(i).vision = [];
    glmtrial(i).audition = [];
    if Z(i,2)==1 % attend to vision rule
        glmtrial(i).vision = 0;
    elseif Z(i,2)==2 % attend to audition rule
        glmtrial(i).audition = 0;
    end;
    
    
    glmtrial(i).lowpasscue = [];
    glmtrial(i).highpasscue = [];
    glmtrial(i).uvled = [];
    glmtrial(i).greenled = [];
    
    glmtrial(i).laserOn = [];    
    
    if Z(i,9)== 0 && Z(i,10) == 0
        if Z(i,2)==1 % vision
            glmtrial(i).R1C1 = cueon;
            glmtrial(i).highpasscue = 0;
        elseif Z(i,2)==2 % audition
            glmtrial(i).R2C1 = cueon;
            glmtrial(i).lowpasscue = 0;
        end
    elseif Z(i,9)== 3 && Z(i,10) == 0
        if Z(i,2) == 1 % vision
            glmtrial(i).R1C2 = cueon;
            glmtrial(i).greenled = 0;
        elseif Z(i,2)==2 % audition
            glmtrial(i).R2C2 = cueon;
            glmtrial(i).uvled = 0;
        end
    end;
    
    
    if Z(i,9)== 0 && Z(i,10) == 1
        if Z(i,2)==1 % vision
            glmtrial(i).R1C1 = cueon;
            glmtrial(i).highpasscue = 0;
            glmtrial(i).laserOn = 0;
        elseif Z(i,2)==2 % audition
            glmtrial(i).R2C1 = cueon;
            glmtrial(i).lowpasscue = 0;
             glmtrial(i).laserOn = 0;
        end
    elseif Z(i,9)== 3 && Z(i,10) == 1
        if Z(i,2) == 1 % vision
            glmtrial(i).R1C2 = cueon;
            glmtrial(i).greenled = 0;
             glmtrial(i).laserOn = 0;
        elseif Z(i,2)==2 % audition
            glmtrial(i).R2C2 = cueon;
            glmtrial(i).uvled = 0;
             glmtrial(i).laserOn = 0;
        end
    end
    
    glmtrial(i).reward = Z(i,5);
    
    for n = 1:numPFC
        temp_spiketimes = pfc_spiketimes(n).SpikeTimes{i}*1000 - trial_start_ms;
        glmtrial(i).(['PFCUnit' num2str(n)]) = temp_spiketimes(temp_spiketimes>=0 & temp_spiketimes<glmtrial(i).duration);
    end
    
    for n = 1:numMD
        temp_spiketimes = md_spiketimes(n).SpikeTimes{i}*1000 - trial_start_ms;
        glmtrial(i).(['MDUnit' num2str(n)]) = temp_spiketimes(temp_spiketimes>=0 & temp_spiketimes<glmtrial(i).duration);
    end
    
end;

%% ==== Compute CMI

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

% CMIPFC = ones( 1,length(keptPFC) );
% CMIMD  = ones( 1,length(keptMD) );


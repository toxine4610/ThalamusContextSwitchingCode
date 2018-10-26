

clear all;

% Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/BPM1_3/ContextSwitchHalo/';
% ExptDates = {'2018-02-11','2018-02-12','2018-02-13', '2018-02-14', '2018-02-16'};
%  ExptDates = {'2018-02-17','2018-02-18','2018-02-20','2018-02-22','2018-02-23','2018-02-24'};

Datapath = '/Users/rvrikhye/Dropbox (Personal)/Rajeev/ForContextSwitchProject/DoubleCueDatabase/BPM1_3/ContextSwitchHalo/UniLatMD-Low/';
ExptDates = {'2018-03-13','2018-03-14','2018-03-15','2018-03-16','2018-03-22'};

% Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/SOMCre/DualModalityCue/';
% ExptDates = {'2017-12-15','2017-12-16','2017-12-17','2017-12-18','2017-12-20','2017-12-22','2017-12-23','2017-12-24','2017-12-26'};

% ,'2018-02-07' , ...
%     '2018-02-10', '2018-02-14', '2018-02-18', '2018-02-20', '2018-02-22', '2018-02-24'


firingRatesAverage = [];

for d = 1:length(ExptDates)
    
 
    clear Spfc Spfc Spfc_mix Spfc_mix Z_C1 D foo
    clear R1_sound_R R1_sound_L R2_sound_R R2_sound_L
    clear R1_nosound_R R1_nosound_L R2_nosound_R R2_nosound_L
    
    load( [Datapath ExptDates{d} '/BPM13_uniMDHalo_Session.mat'] );
%     Z_C1(:,10) = 0;
    
    indCorrVis = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 1 );
    indCorrAud = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 2 );
    
    ZVis = Z_C1(indCorrVis, :);
    ZAud = Z_C1(indCorrAud, :);
    
%     ZVis = ZVis( , :);
%     ZAud = ZAud( find(ZAud(:,10)==1), :);
    
    indVisNoSound_RightChoice = find( ZVis(:,9) == 3 & ZVis(:,10)==1);
    indVisNoSound_LeftChoice = find( ZVis(:,9) == 3 & ZVis(:,10)==0);
    
    indAudNoSound_RightChoice = find( ZAud(:, 9) == 3 & ZAud(:,10)==1);
    indAudNoSound_LeftChoice = find( ZAud(:, 9) == 3 & ZAud(:,10)==0);
    
    indVisSound_RightChoice  = find( ZVis(:,9) == 0 & ZVis(:,10)== 1);
    indVisSound_LeftChoice  = find( ZVis(:,9) == 0 & ZVis(:,10)== 0);
    
    indAudSound_RightChoice   = find( ZAud(:, 9) == 0 & ZAud(:,10)== 1);
    indAudSound_LeftChoice   = find( ZAud(:, 9) == 0 & ZAud(:,10)== 0);
    
    %%
    
    range = [ -0.2, 0.7 ];
    bin = 0.0010;
    filtWidth = 0.08;
    
    
    for i = 1:numel(Spfc)
        
        Spfc_mix(i).SpikeTimes_R1_sound_R = Spfc(i).SpikeTimes_R1C1(indVisSound_RightChoice);
        [R1_sound_R(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R1_sound_R, range , bin, filtWidth);
        
        Spfc_mix(i).SpikeTimes_R1_sound_L = Spfc(i).SpikeTimes_R1C1(indVisSound_LeftChoice);
        [R1_sound_L(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R1_sound_L, range , bin, filtWidth);
        
        Spfc_mix(i).SpikeTimes_R2_sound_R = Spfc(i).SpikeTimes_R2C1(indAudSound_RightChoice);
        [R2_sound_R(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R2_sound_R, range , bin, filtWidth);
        
        Spfc_mix(i).SpikeTimes_R2_sound_L = Spfc(i).SpikeTimes_R2C1(indAudSound_LeftChoice);
        [R2_sound_L(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R2_sound_L, range , bin, filtWidth);
        
        Spfc_mix(i).SpikeTimes_R1_nosound_R = Spfc(i).SpikeTimes_R1C1(indVisNoSound_RightChoice);
        [R1_nosound_R(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R1_nosound_R, range , bin, filtWidth);
        
        Spfc_mix(i).SpikeTimes_R1_nosound_L = Spfc(i).SpikeTimes_R1C1(indVisNoSound_LeftChoice);
        [R1_nosound_L(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R1_nosound_L, range , bin, filtWidth);
        
        Spfc_mix(i).SpikeTimes_R2_nosound_R = Spfc(i).SpikeTimes_R2C1(indAudNoSound_RightChoice);
        [R2_nosound_R(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R2_nosound_R, range , bin, filtWidth);
        
        Spfc_mix(i).SpikeTimes_R2_nosound_L = Spfc(i).SpikeTimes_R2C1(indAudNoSound_LeftChoice);
        [R2_nosound_L(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R2_nosound_L, range , bin, filtWidth);
        
    end;
    
    
    %%
    foo(:, 1, 1, :) = R1_nosound_L;
    foo(:, 1, 2, :) = R1_nosound_R;
    
    foo(:, 2, 1, :) = R2_nosound_L;
    foo(:, 2, 2, :) = R2_nosound_R;
    
    foo(:, 3, 1, :) = R1_sound_L;
    foo(:, 3, 2, :) = R1_sound_R;
    
    foo(:, 4, 1, :) = R2_sound_L;
    foo(:, 4, 2, :) = R2_sound_R;
    
    firingRatesAverage = cat(1, firingRatesAverage, foo);
    
end;

%%
combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

timeEvents = [0, 0.6];

%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% minimal plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);


%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
   'combinedParams', combinedParams);

%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);


%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

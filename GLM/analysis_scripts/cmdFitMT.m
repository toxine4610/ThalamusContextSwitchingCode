% Fit GLM to MT units
% This code calls the fitting code and analysis code for each LIP unit in
% the dataset. The fits are regularized using ridge-regression. This code
% uses a fixed ridge parameter (set below)
% All analyses are handled by the mtlipglm class. The design/fitting of
% each GLM is handled by the glmspike class, which is an instance of the
% neuroGLM class available at https://github.com/jcbyts/neuroGLM
dataPath = getpref('mtlipglm', 'dataPath');

% --- Find experiments that have multiple MT units
experiments=getExperimentsAnd({'MT', 'simultaneous'});
regenerateModelComparisonFiles = true; % regenerate model comparison files from fits
ridgeParameter = 0.1; % regularize fits with ridge-regression (L2 penatly)
refitModel     = false; % if true, fits will be overwritten

% Parfor if possible, change to for-loop if you do not have the parallel
% toolbox
parfor kExperiment=1:numel(experiments)
    
    % --- setup analysis for each session
    exname=experiments{kExperiment};
    mstruct=mtlipglm(exname, dataPath, 'overwrite', refitModel, 'nFolds', 5);
    mstruct.overwrite =refitModel;
    mstruct.binSize = 10; % Don't change. Analysis code assumes binSize=10
    mstruct.modelDir = fullfile(dataPath,'main_fits');
    % buildTrialStruct builds a struct array of each trial's data in the
    % format that is used by neuroGLM. This is way faster than the LIP fits
    % because there is no need to simulate MT input ('MTsim', false)
    mstruct.buildTrialStruct('IncludeContrast', true, 'MotionEpoch', true, 'MTsim', false);
    
    for kNeuron=1:numel(mstruct.neurons)
        if ~mstruct.neurons(kNeuron).isMT % skip LIP
            continue
        end
        fname=fullfile(mstruct.modelDir, [mstruct.neurons(kNeuron).getName 'modelComparison.mat']);
        
        if exist(fname, 'file')&&~mstruct.overwrite&&~regenerateModelComparisonFiles
            fprintf('\n\n\nSkipping [%s] \n\n\n\n\n', fname)
            continue
        else
            % --- Fitting happens here
            % Specify which models to fit:
            % 1 - Stimulus-to-LIP
            % 2 - Stimulus-to-LIP w/ history filter
            % 3 - Inter-Area Coupling
            % 4 - Intra-Area Coupling
            modelIx = 1:4;
            g=fitAllModels(mstruct,kNeuron,true,modelIx,'instantaneousCoupling', true, 'rho', ridgeParameter);
            P=mstruct.modelComparisonToo(g);
            parsave(fname, P,'-v7')
        end
        
    end
end

pushMessage('Done running cmdFitModels.m')

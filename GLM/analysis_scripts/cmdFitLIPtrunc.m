% Fit GLM to LIP units
% This code calls the fitting code and analysis code for each LIP unit in
% the dataset. The fits are regularized using ridge-regression. This code
% uses a fixed ridge parameter (set below)
% All analyses are handled by the mtlipglm class. The design/fitting of
% each GLM is handled by the glmspike class, which is an instance of the
% neuroGLM class available at https://github.com/jcbyts/neuroGLM
dataPath = getpref('mtlipglm', 'dataPath');

% --- Find experiments that have multiple LIP units
experiments=getExperimentsAnd({'LIP', 'simultaneous'});
regenerateModelComparisonFiles = true; % regenerate model comparison files from fits
ridgeParameter = 0.1; % regularize fits with ridge-regression (L2 penatly)
refitModel     = false; % if true, fits will be overwritten
truncations    = 100:100:3e3; % time in ms preceding saccade
% Parfor if possible, change to for-loop if you do not have the parallel
% toolbox
parfor kExperiment=1:numel(experiments)
    
    
    % --- setup analysis for each session
    exname=experiments{kExperiment};
    mstruct=mtlipglm(exname, dataPath, 'overwrite', refitModel, 'nFolds', 5);
    mstruct.overwrite =refitModel;
    mstruct.binSize = 10; % Don't change. Analysis code assumes binSize=10
    mstruct.modelDir = fullfile(dataPath,'lip_trunc_fits');
    % buildTrialStruct builds a struct array of each trial's data in the
    % format that is used by neuroGLM. This is very slow when the MTsim
    % flag is turned on because it has to simulate from the MT model on
    % each trial
    mstruct.buildTrialStruct('IncludeContrast', true, 'MotionEpoch', true, 'MTsim', false);
    
    for k=1:numel(mstruct.trial)
        mstruct.trial(k).gaborContrast=mstruct.trial(k).contrasts;
    end
    
    isLIP=find([mstruct.neurons(:).isLIP]);
    
    for kNeuron=isLIP(:)'
        
        % filename for saved analyses
        fname=fullfile(mstruct.modelDir, [mstruct.neurons(kNeuron).getName 'modelComparison.mat']);
        
        if exist(fname, 'file')&&~mstruct.overwrite&&~regenerateModelComparisonFiles
            fprintf('\n\n\nSkipping [%s] \n\n\n\n\n', fname)
            continue
        else
            g=mstruct.fitLIPtruncation(kNeuron, truncations, rho)
            S=mstruct.modelComparisonToo(g);
            parsave(fname, S,'-v7')
        end
        
    end
end

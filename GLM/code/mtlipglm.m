classdef mtlipglm < handle
    % MTLIPGLM class for handling experimental data from MT and LIP with the
    % gabor pulse task
    %
    %  mstruct=mtlipglm();
    %     m.overwrite=false;
    %     m.nFolds=5;
    %     m.regressionMode='RIDGEFIXED';
    %     m.binSize=10;
    %     m.modelDir=modelDir;
    %     m.getExperiment(exname)
    %     m.buildTrialStruct
    
    
    properties
        exname@char             % experiment name
        stim@rcgreader          % stimulus object
        neurons@neuro           % neuron objects
        lfp                     % local field potential
        directory               % path to where data lives
        trial@struct            % struct array of trial data
        trialParam@struct       % experiment parameters
        stimValidTrials@double  % index for trials in the stimulus object
        nFolds=5                % number of CV folds for evaluating model
        overwrite@logical=false % refit the models and overwrite files
        modelDir                % name of directory to save/load fits
        binSize@double=10;      % bin size in ms
    end
    
    methods
        %% Constructor
        function obj = mtlipglm(exname, dataPath, varargin)
            
            p=inputParser();
            p.addOptional('overwrite', obj.overwrite, @islogical)
            p.addOptional('nFolds', obj.nFolds, @isnumeric)
            p.parse(varargin{:});
            
            obj.overwrite   = p.Results.overwrite;
            obj.nFolds      = p.Results.nFolds;
            obj.directory   = dataPath;
            
            obj.getExperiment(exname);
        end % constructor
        
        % --- Get Neuron/Stimulus Data
        function getExperiment(obj,exname)
            % getExperiment(exname)
            obj.exname = exname;
            obj.stim=rcgreader(obj.exname, obj.directory);
            obj.neurons=getNeurons(obj.exname,obj.directory);
        end
        
        % --- Build Trial Struct for NeuroGLM
        function buildTrialStruct(obj, varargin)
            % BUILDTRIALSTRUCT
            % Builds struct array of trial data for NeuroGLM
            % Optional Arguments include:
            % IncludeEyepos@logical   - include eye posisiont, pupil,
            %                           saccade start/end
            % IncludeLFP@logical      - include local field potential
            %                           (requires LFP data files in path)
            % IncludeContrast@logical - use phase of gabors to compute
            %                           contrast (opposed to just using a
            %                           boxcar for the stimulus value)
            % MotionEpoch@logical     - truncate trial to valid fixation
            % MTsim@logical           - simulate MT rates as stimulus
            %                           (needed for MT-to-LIP model)
            
            
            for k=1:numel(obj.neurons)
                nw(k)=obj.neurons(k).getStruct(); %#ok<AGROW>
            end
            
            st=getStim(obj.exname, obj.directory);
            
            fprintf('Building Trial Structure [%s]\n', st.exname)
            
            p = inputParser();
            p.addOptional('IncludeContrast', false)
            p.addOptional('IncludeLFP', false)
            p.addOptional('IncludeEyepos', false)
            p.addOptional('MotionEpoch', false)
            p.addOptional('MTsim',false)
            p.parse(varargin{:});
            
            if isempty(obj.lfp)
                includeLFP = false;
            else
                includeLFP = p.Results.IncludeLFP;
            end
            
            binres      = 1e-3; % discretize everything at 1ms resolusion
            binfun      = @(t) (t==0) + ceil(t/binres);
            nNeurons    = numel(nw);
            
            
            if includeLFP
                fprintf('Found %d nw across %d channels of LFP\n', nNeurons, numel(obj.lfp.info.channels))
            end
            
            
            assert(binres==1e-3, 'Binning has only been tested at 1ms bins.')
            
            %--------------------------------------------------------------
            % To simulate MT responses, use filters recovered from GLM fits
            % to the MT neurons in our population
            if p.Results.MTsim
                load(fullfile(obj.directory, 'MTkernels.mat'))
                
                % upsample filters to get 1ms resolution
                kContr1 = upsampleFilters(kCons, obj.binSize);
                kDir1   = upsampleFilters(kDirs, obj.binSize);
                
                % simulate MT responses for each monkey
                if strcmp(obj.exname(1), 'p')
                    kContr1 = kContr1(:,isPat);
                    kDir1   = kDir1(:,isPat);
                    b0      = b0(isPat);
                else
                    kContr1 = kContr1(:,~isPat);
                    kDir1   = kDir1(:,~isPat);
                    b0      = b0(~isPat);
                end
            end
            
            % Find trials that are good for all neurons
            ndx=false(st.nTrials, nNeurons);
            for kNeuron = 1:nNeurons
                ndx(nw(kNeuron).trialIndex,kNeuron) = true;
            end
            
            % remove neurons with fewer than 100 trials
            remix = sum(ndx,1) < 100;
            if any(remix)
                error('Some neurons are removed')
            end
            ndx(:,remix) = [];
            
            %--------------------------------------------------------------
            % Setup variables to include
            % find valid trials (good behavior and have stable neural recordings)
            trialStartsPlexonTime = [st.timing(:).plxstart]';
            validTrials = find(all(ndx,2) & st.goodtrial(:) & ~isnan(trialStartsPlexonTime));
            nTrials = numel(validTrials);
            
            timingFields = {'duration','motionon', 'motionoff', 'fpon', 'fpentered', 'gosignal', 'targson', 'saccade', 'pulseon', 'reward'};
            
            valueFields = {'choice', 'pulses', 'trialNumber', 'trialId', 'isRevco'};
            
            if p.Results.IncludeEyepos
                valueFields = {'eyepos',  'saccadeStarts', 'saccadeEnds', valueFields{:}}; %#ok<CCAT>
                eyepos=st.eyepos;
                if abs(mode(diff(eyepos{1}(:,1)))-.001)>.001
                    eyepos=upsampleEyepos(eyepos, 1/binres);
                end
            end
            
            if includeLFP
                valueFields = {'LFP',  valueFields{:}}; %#ok<CCAT>
                trialStartLFPSamples = plx.convertTimeToSamples(trialStartsPlexonTime, obj.lfp.info.adfreq, obj.lfp.info.timestamps, obj.lfp.info.fragsamples);
            end
            
            % --- neuron names
            neuronFields  = cell(1,nNeurons);
            neuronTrialDC = cell(1,nNeurons);
            for kNeuron = 1:nNeurons
                %                 neuronFields{kNeuron} = nw(kNeuron).getName(false);
                neuronFields{kNeuron} = sprintf('%sneuron%02.0fch%02.0f', nw(kNeuron).brainArea, nw(kNeuron).id, nw(kNeuron).channel);
                neuronTrialDC{kNeuron} = ['FR_' neuronFields{kNeuron}];
            end
            
            % --- Initialize struct array
            fields = [timingFields valueFields neuronFields neuronTrialDC];
            fields = [fields; cell(1, numel(fields))];
            
            obj.trial = repmat(struct(fields{:}), nTrials, 1);
            
            % --- Check which target correlates with positive pulses
            chocorr=corr(sum(sum(st.pulses(st.goodtrial,:,:),2),3), st.targchosen(st.goodtrial));
            if chocorr<0
                tpositive=1;
            else
                tpositive=2;
            end
            
            %--------------------------------------------------------------
            % --- Loop over trials and build up values
            for kTrial = 1:nTrials
                t = validTrials(kTrial);
                
                % --- timing variables - convert to ms bins
                obj.trial(kTrial).duration  = binfun(st.timing(t).duration);
                obj.trial(kTrial).motionon  = binfun(st.timing(t).motionon);
                obj.trial(kTrial).motionoff = binfun(st.timing(t).motionoff);
                obj.trial(kTrial).pulseon   = binfun(st.timing(t).pulses);
                obj.trial(kTrial).pulseoff  = [binfun(st.timing(t).pulses(2:end)); binfun(st.timing(t).motionoff)];
                obj.trial(kTrial).saccade   = binfun(st.timing(t).choice);
                obj.trial(kTrial).fpon      = binfun(st.timing(t).fpon);
                obj.trial(kTrial).fpentered = binfun(st.timing(t).fpentered);
                obj.trial(kTrial).gosignal  = binfun(st.timing(t).fpoff);
                obj.trial(kTrial).isRevco   = st.dirprob(t)==0;
                if isempty(obj.trial(kTrial).gosignal)
                    obj.trial(kTrial).gosignal=nan;
                end
                obj.trial(kTrial).targson   = binfun(st.timing(t).targson);
                obj.trial(kTrial).reward    = binfun(st.timing(t).reward);
                
                % --- in some sessions targets were turned off at reward
                if isfield(st.timing(t), 'targsoff') && ~isnan(st.timing(t).targsoff)
                    obj.trial(kTrial).targsoff  = binfun(st.timing(t).targsoff);
                else
                    obj.trial(kTrial).targsoff  = binfun(st.timing(t).reward);
                end
                
                % --- Value fields
                obj.trial(kTrial).choice  = st.targchosen(t)==tpositive;
                obj.trial(kTrial).pulses  = mean(st.pulses(t,:,:),3);
                obj.trial(kTrial).trialNumber=st.trialnumber(t);
                obj.trial(kTrial).trialId=st.trialId(t);
                
                % --- Motion stimulus
                frameInds = obj.trial(kTrial).motionon:obj.trial(kTrial).motionoff;
                
                % boxcar for when motion is on
                obj.trial(kTrial).motion = sparse(frameInds, ones(1, numel(frameInds)), ones(1, numel(frameInds)), obj.trial(kTrial).duration, 1);
                
                % include parameters derived from individual gabor values
                if p.Results.IncludeContrast
                    
                    % get contrast on each frame and change in phase
                    [contrast, ~, dphi] = getSingleGabor(squeeze(st.maskphase(t,:,:)), squeeze(st.basephase(t,:,:)));
                    
                    % average phase change on each frame
                    vals = mean(squeeze(contrast).^2 .* squeeze(dphi),2);
                    iis  = round(linspace(1, numel(vals), numel(frameInds)));
                    vals = vals(iis);
                    obj.trial(kTrial).direction=sparse(frameInds, ones(1, numel(frameInds)), vals, obj.trial(kTrial).duration, 1);
                    
                    % average "contrast" of each gabor
                    contrast=squeeze(contrast).^2;
                    dphi=squeeze(dphi);
                    nGabors=size(contrast,2);
                    iiv=bsxfun(@plus, iis', (0:nGabors-1)*size(contrast,1));
                    iix=repmat(frameInds(:), 1, nGabors);
                    iiy=ones(numel(frameInds),1)*(1:nGabors);
                    obj.trial(kTrial).contrasts=sparse(iix(:), iiy(:), contrast(iiv(:)), obj.trial(kTrial).duration, nGabors);
                    obj.trial(kTrial).deltaphase=sparse(iix(:), iiy(:), dphi(iiv(:)), obj.trial(kTrial).duration, nGabors);
                    spavcon=mean(contrast,2); % spatially averaged and thresholded contrast
                    
                    obj.trial(kTrial).contrast=sparse(frameInds(:), ones(numel(frameInds),1), spavcon(iis), obj.trial(kTrial).duration, 1);
                    
                    % simulate from MT population
                    if p.Results.MTsim
                        
                        % --- Simulate from full population
                        tic % this is slow
                        Xdir=basisFactory.boxcarStim(obj.trial(kTrial).pulseon, obj.trial(kTrial).pulseoff, obj.trial(kTrial).duration, obj.trial(kTrial).pulses);
                        Xdir = full(Xdir);
                        
                        if contrastIsBoxcar % if using boxcar
                            contrast = obj.trial(kTrial).motion;
                        else % if using full contrast
                            contrast = full(obj.trial(kTrial).contrast)/10;
                            contrast = contrast(:);
                        end
                        
                        pulse    = Xdir(:);
                        
                        % simulate population of MT responses
                        C = convmtx(contrast,size(kContr1,1));
                        D = convmtx(pulse, size(kDir1,1));
                        MTp=exp(bsxfun(@plus, b0, C*kContr1 + D*kDir1));
                        MTa=exp(bsxfun(@plus, b0, C*kContr1 + D*-kDir1));
                        
                        MTp = mean(MTp,2);
                        MTa = mean(MTa,2);
                        
                        MTp = MTp-MTp(1); % subtract baseline FR
                        MTa = MTa-MTa(1);
                        
                        obj.trial(kTrial).MTpref=MTp(1:obj.trial(kTrial).duration);
                        obj.trial(kTrial).MTanti=MTa(1:obj.trial(kTrial).duration);
                        
                        et = toc;
                        fprintf('Trial: %d took %d ms\n', kTrial, et*1e3)
                        
                        
                    end % simulate MT
                end % include Contrast
                
                
                % -- Trial spike times
                for kNeuron = 1:nNeurons
                    tt = (nw(kNeuron).spikeTimes(nw(kNeuron).spikeTimes > trialStartsPlexonTime(t) & nw(kNeuron).spikeTimes < (trialStartsPlexonTime(t) + st.timing(t).duration)) - trialStartsPlexonTime(t));
                    obj.trial(kTrial).(neuronFields{kNeuron}) = binfun(tt);
                    
                    % --- pre-trial spike rate
                    obj.trial(kTrial).(neuronTrialDC{kNeuron}) = sum( (nw(kNeuron).spikeTimes > trialStartsPlexonTime(t) - 1) & (nw(kNeuron).spikeTimes < trialStartsPlexonTime(t)));
                end
                
                % --- Include Eye Data
                if p.Results.IncludeEyepos
                    % eyepos
                    trialEyepos=eyepos{t};
                    eyendx = (trialEyepos(:,1) >= 0 & trialEyepos(:,1) <= st.timing(t).duration);
                    [result, smoothTrace] = saccadeDetector(trialEyepos(eyendx,1)', trialEyepos(eyendx,2:3)');
                    obj.trial(kTrial).eyepos = smoothTrace';
                    obj.trial(kTrial).saccadeStarts = ceil((result(1,:)-trialEyepos(find(eyendx,1),1))*1e3);
                    obj.trial(kTrial).saccadeEnds   = ceil((result(2,:)-trialEyepos(find(eyendx,1),1))*1e3);
                    
                    if size(trialEyepos,2) == 4
                        obj.trial(kTrial).pupil  = trialEyepos(eyendx, 4);
                    end
                end
                
                % --- Local Field Potentials
                if includeLFP
                    idx = trialStartLFPSamples(t) + (0:(ceil(st.timing(t).duration*obj.lfp.info.adfreq)-1));
                    obj.trial(kTrial).LFP = obj.lfp.data(idx, :);
                end
            end
            
            % --- Truncate trial to epoch after fixation is obtained (to
            % eliminate as many uncontrolled eye movements as possible)
            if p.Results.MotionEpoch
                fields=fieldnames(obj.trial);
                
                % only keep timing fields
                fields = setdiff(fields, valueFields); % ignore valiue fields
                fields = setdiff(fields, neuronTrialDC); % ignore trial DC firing rates
                
                for kTrial=1:numel(obj.trial)
                    
                    mo=obj.trial(kTrial).fpentered;
                    go=obj.trial(kTrial).saccade;
                    
                    % count spike rate before trial start
                    for kNeuron = 1:numel(neuronFields)
                        obj.trial(kTrial).(neuronTrialDC{kNeuron}) = sum(obj.trial(kTrial).(neuronFields{kNeuron})<mo)/mo*1e3;
                    end
                    
                    d=go-mo;
                    for kField=1:numel(fields)
                        if strcmp(fields{kField}, 'duration') || strcmp(fields{kField}, 'choice') || strcmp(fields{kField}, 'pulses')
                            continue
                        end
                        
                        % check if continuous
                        if any(size(obj.trial(kTrial).(fields{kField}))==obj.trial(kTrial).duration)
                            obj.trial(kTrial).(fields{kField})=obj.trial(kTrial).(fields{kField})(mo:go-1,:);
                            continue
                        end
                        
                        % adjust timing
                        t=obj.trial(kTrial).(fields{kField})-mo;
                        t(t>d | t<0)=[];
                        obj.trial(kTrial).(fields{kField})=t;
                        
                    end % timing fields loop
                    
                    if isempty(obj.trial(kTrial).gosignal)
                        obj.trial(kTrial).gosignal = d;
                    end
                    
                    obj.trial(kTrial).duration=d;
                    obj.trial(kTrial).motiononly=true;
                end % trial loop
                
                
            end % if motion epoch
            
            
            % --- build params struct
            pfields = {'nPulses', 'luminanceBackground', 'luminanceMinMuMax', 'gaborXY', 'sf', 'theta', 'tf', 'nFrames', 'nGabors', 'frate', 'ppd'};
            args = [pfields; cell(1, numel(pfields))];
            obj.trialParam = struct(args{:});
            for k = 1:numel(pfields)
                obj.trialParam.(pfields{k}) = st.(pfields{k});
            end
            obj.trialParam.pulsedur = st.nFrames/st.nPulses/st.frate*1e3;
            obj.trialParam.stimValidTrials = validTrials;
            obj.trialParam.trialCnt=st.trialCnt;
            obj.stimValidTrials=validTrials;
        end
        
        % --- fit all the models
        function g=fitAllModels(obj, kNeuron, saveIt, modelIndex, varargin)
            % FITALLMODELS fits all the models described in Yates et al., 2017
            % glm = fitAllModels(obj, kNeuron, saveIt, modelIndex, varargin)
            
            
            % parse arguments to find hyperparameter(s)
            args = varargin;
            rhoarg = find(cellfun(@(x) strcmp(x, 'rho'), args));
            if ~isempty(rhoarg)
                rho=args{rhoarg+1};
                args([rhoarg rhoarg+1])=[];
            else
                rho = 0.1; % default to small amount of regularization
            end
            
            if ~exist('saveIt', 'var') || isempty(saveIt)
                saveIt=true;
            end
            
            % --- LIP has more sluggish responses than MT so the bases must
            % span a wider time window. <intParams> sets the time window to
            % analyze for certain regressors
            % intParams = [motion, history, couplingInter, couplingIntra]
            % (all intParams are in miliseconds)
            if obj.neurons(kNeuron).isMT
                intParams=[400, 200, 500, 500];
            else
                intParams=[800, 200, 500, 500];
            end
            
            fprintf('********************************************************\n')
            fprintf('********************************************************\n')
            fprintf('setting up fit [%s]\n', obj.neurons(kNeuron).getName(0))
            
            modelNames={'Poisson', 'Uncoupled', 'InterArea', 'IntraArea', 'MTsim', 'MTsimChoice'};
            if ~exist('modelIndex', 'var')||isempty(modelIndex)
                modelIndex=1:numel(modelNames);
            end
            
            nModels=numel(modelIndex);
            mparams={[1 0 0 0],[1 1 0 0], [1 1 0 1], [1 1 1 0], [1 0 0 0], [1 0 0 0]};
            
            % -- loop over models
            for kModel=1:nModels
                idModel=modelIndex(kModel);
                
                % --- Load fit if the file exists
                if exist(fullfile(obj.modelDir, [obj.neurons(kNeuron).getName modelNames{idModel} '.mat']), 'file') && ~obj.overwrite
                    fprintf('********************************************************\n')
                    fprintf('********************************************************\n')
                    fprintf('File Exists. Loading and skipping [%s]\n', modelNames{idModel})
                    fprintf('********************************************************\n')
                    fprintf('********************************************************\n')
                    tmp=load(fullfile(obj.modelDir, [obj.neurons(kNeuron).getName modelNames{idModel}]));
                    % this is to deal with disappearing function handle
                    % bugs in certain versions of matlab
                    g(kModel)=struct2funh(funh2struct(tmp.obj,false),false);
                    continue
                end
                
                % --- Fit the model
                g(kModel)=glmspike(modelNames{idModel});
                
                % set parameters that govern the design matrix
                g(kModel).setDesignOptions('binSize', obj.binSize, ...
                    'boxcarPulse', true, args{:});
                
                % Turn this on to use a boxcar to represent the stimulus
                % contrast. Otherwise, this will use the spatially averaged
                % contrast of the gabors. If false, IncludeContrst must be
                % true during buildTrialStruct()
                g(kModel).setDesignOptions('boxcarContrast', false);
                
                % If fitting a coupling model or history model, make sure
                % that any residual error from the onset transient is not
                % absorbed by the coupling filters
                if strncmp(modelNames{idModel}, 'In',2) || strncmp(modelNames{idModel}, 'Un',2)
                    g(kModel).setDesignOptions('includeMotionOnset', true);
                end
                
                % Turn on for the MT-to-LIP model
                if any(strfind(modelNames{idModel}, 'MTsim'))
                    g(kModel).setDesignOptions('MTsim', 1)
                end
                
                % Turn on for full choice term instead of truncated choice
                % term
                if any(strfind(modelNames{idModel}, 'Choice'))
                    g(kModel).setDesignOptions('includeChoice', true);
                end
                
                % Final setting of parameters
                g(kModel).param=obj.trialParam;
                g(kModel).fitNeuron=obj.neurons(kNeuron).getName(0);
                
                % turn off unused stimulus kernels
                params=intParams.*mparams{idModel};
                
                % create temporal basis functions
                g(kModel).buildDesignBoxcar(obj.trial, params(1), params(2), params(3), params(4));
                
                % build design matrix
                g(kModel).compileDesignMatrix(obj.trial, 1:numel(obj.trial))
                
                % add a column of ones to capture the baseline firing rate
                g(kModel).addBiasColumn('right')
                
                % set method for learning hyperparameters (disabled in
                % release code)
                g(kModel).model.regressionMode='RIDGEFIXED';
                
                
                g(kModel).fitCV(obj.trial, obj.nFolds, rho);
                
                if saveIt
                    g(kModel).save(obj.modelDir, obj.exname)
                    g(kModel)=struct2funh(g(kModel), false);
                end
                
            end % Loop over models
            
        end % fitAllModels
        
        % --- Fit the LIP model with varying truncations of choice term
        function g=fitLIPtruncation(obj, kNeuron, truncations, rho)
            % FITLIPTRUNCATION fits the stim-to-LIP model with varying
            % amounts of pre-saccadic term
            % g=fitLIPtruncation(obj, kNeuron, truncations, rho)
            % Inputs:
            %   obj@mtlipglm    - can be called m.fitLIPtruncation(...)
            %   kNeuron@double  - neuron number
            %   truncations@double - array of times pre-saccade (in ms)
            %   rho@double      - hyper-parameter
            
            % --- loop over truncations
            for kModel=1:numel(truncations)
                modelName=sprintf('sac%d', truncations(kModel));
                fname=fullfile(obj.modelDir, [obj.neurons(kNeuron).getName sprintf('sac%d', truncations(kModel)) '.mat']);
                if exist(fname, 'file') && ~obj.overwrite
                    fprintf('********************************************************\n')
                    fprintf('********************************************************\n')
                    fprintf('File Exists. Loading and skipping [%s]\n', modelName)
                    fprintf('********************************************************\n')
                    fprintf('********************************************************\n')
                    tmp=load(fname);
                    g(kModel)=struct2funh(funh2struct(tmp.obj,false),false);
                    continue
                end
                
                g(kModel)=glmspike(modelName);
                g(kModel).setDesignOptions('binSize', obj.binSize, 'boxcarPulse', true);
                g(kModel).setDesignOptions('includeChoice', -truncations(kModel));
                
                g(kModel).param=obj.trialParam;
                g(kModel).fitNeuron=obj.neurons(kNeuron).getName(0);
                params=[800, 0, 0,0];
                g(kModel).buildDesignBoxcar(obj.trial, params(1), params(2), params(3), params(4));
                g(kModel).compileDesignMatrix(obj.trial, 1:numel(obj.trial))
                g(kModel).addBiasColumn('right')
                g(kModel).model.regressionMode='RIDGEFIXED';
                

                g(kModel).fitCV(obj.trial, obj.nFolds, rho);

                
                g(kModel).save(obj.modelDir, obj.stim.exname);
                g(kModel)=struct2funh(g(kModel), false);
            end
            
        end
        
        % --- Check if LIP units need a direction term
        % Is the direction kernel necessary? What percent of the population
        % have their fits improved with the inclusion of a direction kernel
        function g=fitLIPStim(obj, kNeuron, rho)
            
            % there are two models: 1 with pulses. 1 without
            modelNames{1}='NoPulse';
            modelNames{2}='Pulse';
            modelOpts{1} = {'excludeMotion', true, 'includeContrast', true, 'includeChoice', -2.5e3};
            modelOpts{2} = {'excludeMotion', false, 'includeContrast', true, 'includeChoice', -2.5e3};
            
            
            for kModel=1:numel(modelOpts)
                modelName=modelNames{kModel};
                fname=fullfile(obj.modelDir, [obj.neurons(kNeuron).getName sprintf('directionDependncy_%s', modelNames{kModel}) '.mat']);
                if exist(fname, 'file') && ~obj.overwrite
                    fprintf('********************************************************\n')
                    fprintf('********************************************************\n')
                    fprintf('File Exists. Loading and skipping [%s]\n', modelName)
                    fprintf('********************************************************\n')
                    fprintf('********************************************************\n')
                    tmp=load(fname);
                    g(kModel)=struct2funh(funh2struct(tmp.obj,false),false);
                    continue
                end
                g(kModel)=glmspike(modelName);
                g(kModel).setDesignOptions('binSize', obj.binSize, 'boxcarPulse', true);
                g(kModel).setDesignOptions(modelOpts{kModel}{:});
                
                g(kModel).param=obj.trialParam;
                g(kModel).fitNeuron=obj.neurons(kNeuron).getName(0);
                params=[800, 0, 0,0];
                g(kModel).buildDesignBoxcar(obj.trial, params(1), params(2), params(3), params(4));
                g(kModel).compileDesignMatrix(obj.trial, 1:numel(obj.trial))
                g(kModel).addBiasColumn('right')
                g(kModel).model.regressionMode=obj.regressionMode;
                

                g(kModel).fitCV(obj.trial, obj.nFolds, rho);

                
                g(kModel).save(obj.modelDir, obj.stim.exname);
                g(kModel)=struct2funh(g(kModel), false);
            end
            
        end
        
        % --- model comparison
        function S=modelComparison(obj,g)
            % S=modelComparison(obj,g)
            % takes in glmspike objects and evaluates the goodness-of-fit,
            % and other model comparison metrics
            
            assert(obj.binSize==10, 'bin size is assumed tp be 10ms in this code')
            
            win=[-100 300]; % assuming 10ms bins!
            
            nGabors = obj.trialParam.nGabors; % normalizer for PTA
            kModel  = 1; % the Data
            kNeuron = find(arrayfun(@(x) strcmp(x.getName(0), g(kModel).fitNeuron), obj.neurons));
            
            y=g(kModel).getBinnedSpikeTrain(obj.trial, g(kModel).fitNeuron);
            
            sm=10; % 50ms smoothing
            smrate=smooth(full(y),sm)/.01;
            
            dpcell=obj.neurons(kNeuron).dPrimeSigned; %#ok<FNDSB>
            [cohp, psthT,nTrialsPerCoh,Cohs, dpcell]=computeCohPSTH(smrate, g(kModel), obj.trial, win, false, true, dpcell);
            
            % --- Compute the PTA
            ptaWin  = [-20 300];
            ptaBins = 110;
            % basic uncorrected
            [pta, ptaT, ~, dpcell]=computePTA(smrate, g(kModel), obj.trial, 0, ptaWin, ptaBins, dpcell);
            
            % --- Basic structure: filled with the data
            data.y          = y;
            data.psthTime   = psthT;
            data.psthCoh    = cohp;
            data.ptaTime    = ptaT;
            data.ptaRaw     = pta/nGabors;
            data.r2psth     = nan;
            data.llpoisson  = nan;
            data.bps        = nan; % bits / spike
            data.name       = 'data';
            data.wts        = struct;
            data.nTrialsPerCoh = nTrialsPerCoh;
            data.CohEdges   = Cohs;
            data.modelfit   = [];
            data.dprime     = dpcell;
            data.psthFlip   = dpcell;
            data.logliTime  = nan(numel(psthT),obj.nFolds);
            data.smoothingWindow = sm;
            data.designOpts = struct;
            
            % --- the loglikelihood of the data given perfect prediction
            motionon=find(g(kModel).getBinnedSpikeTrain(obj.trial, 'motionon'));
            ybn=binContinuous(y, motionon, win);
            for kFold=1:obj.nFolds
                for kBin=1:numel(psthT)
                    testIndices=g(kModel).modelfit(kFold).testIndices;
                    tmp=ybn(testIndices,kBin);
                    data.logliTime(kBin, kFold)=logliPoisson(tmp(tmp==1), tmp(tmp==1));
                end
            end
            
            % --- What is being compared here (experiment, neuron, models)
            S.exname=obj.stim.exname;
            S.neuron=g(kModel).fitNeuron;
            nModels=numel(g);
            S.model=repmat(data, nModels+1,1);
            
            % --- Loop over models
            for kModel=1:nModels
                
                S.model(kModel+1).name=g(kModel).description;
                S.model(kModel+1).designOpts = g(kModel).designOptions;
                
                % --- correct for the binsize when I used ridge regression
                if strcmp(g(kModel).modelfit(1).regressionMode, 'RIDGE') % RIDGE uses binsize of 1
                    for kFold=1:numel(g(kModel).modelfit)
                        g(kModel).modelfit(kFold).dt=1;
                    end
                end
                
                lambda=g(kModel).cvPredictRate(obj.trial);
                
                smrate=lambda/.01;
                S.model(kModel+1).smoothingWindow = 0;
                
                
                % --- coherence sorted PSTH
                [cohp, psthT, nTrialsPerCond, ~, choFlip]=computeCohPSTH(smrate, g(kModel), obj.trial, win, false, true, data.dprime);
                
                % --- Compute PTA
                % PTA
                [pta, ptaT, ~, dpcell] = computePTA(smrate, g(kModel), obj.trial, 0, ptaWin, ptaBins,data.dprime);

                % --- Update struct with results from this model
                S.model(kModel+1).modelfit  = g(kModel).modelfit;
                S.model(kModel+1).dprime    = dpcell;
                S.model(kModel+1).psthFlip  = choFlip;
                S.model(kModel+1).y         = lambda;
                S.model(kModel+1).psthTime  = psthT;
                S.model(kModel+1).psthCoh   = cohp;
                S.model(kModel+1).ptaTime   = ptaT;
                S.model(kModel+1).ptaRaw    = pta/nGabors;
                
                % --- how good is the model
                tix=psthT>-200& psthT < 1500; % index for evaluating r-squared
                
                mrtrue = S.model(1).psthCoh(tix,nTrialsPerCond>=20); % only take conditions with 20 trials or more
                mrhat  = S.model(kModel+1).psthCoh(tix,nTrialsPerCond>=20);
                
                gix=~isnan(mrtrue) & ~isnan(mrhat); % exclude nans
                
                % calculate rsqaured
                S.model(kModel+1).r2psth    = rsquared(mrtrue(gix), mrhat(gix));
                S.model(kModel+1).llpoisson = logliPoisson(lambda, y);
                S.model(kModel+1).bps       = bitsPerSpike(lambda, y, mean(y)*ones(numel(y),1));
                
                % --- get fitted parameters
                khat = cell2mat(arrayfun(@(x) x.khat, g(kModel).modelfit, 'UniformOutput', false));
                S.model(kModel+1).wts = g(kModel).combineWeights(khat(:,1));
                
                fnames = fieldnames(S.model(kModel+1).wts);
                for kFold=2:size(khat,2)
                    
                    tmpwts=g(kModel).combineWeights(khat(:,kFold));
                    for kField=1:numel(fnames)
                        if strcmp(fnames{kField}, 'bias')
                            S.model(kModel+1).wts.bias(kFold)=tmpwts.bias;
                        else
                            S.model(kModel+1).wts.(fnames{kField}).data=[S.model(kModel+1).wts.(fnames{kField}).data tmpwts.(fnames{kField}).data];
                        end
                    end
                end
                
                % keep track of parameter names
                S.model(kModel+1).covarLabels = {g(kModel).covar(:).label};
                S.model(kModel+1).covarDesc   = {g(kModel).covar(:).desc};
                S.model(kModel+1).khat        = khat;
                
                % --- evaluate the log-likelihood of test data for each bin
                ybn = binContinuous(y, motionon, win); % Data
                lbn = binContinuous(lambda, motionon, win); % Model
                
                % loop over cross-validation test sets
                S.testIndices={g(kModel).modelfit(:).testIndices};
                for kFold=1:obj.nFolds
                    
                    for kBin=1:numel(psthT)
                        testIndices=g(kModel).modelfit(kFold).testIndices;
                        S.model(kModel+1).logliTime(kBin, kFold)=logliPoisson(lbn(testIndices,kBin), ybn(testIndices,kBin));
                    end % time bins
                    
                end % cross validation folds
                
            end % model loop
            
        end % function modelComparisonToo
        
        
        %% plot Maps
        function h=plotMaps(obj)
            h=gcf;
            nNeurons=numel(obj.neurons);
            sx=5;
            sy=5;
            ax=tight_subplot(sy,sx,.002, .002);
            for kNeuron=1:nNeurons
                set(gcf, 'currentaxes', ax(kNeuron))
                obj.neurons(kNeuron).plotMap; axis on
                set(gca, 'XTick', '', 'YTick', '')
                yd=ylim;
                xd=xlim;
                n=[num2str(kNeuron,2) ' (' num2str( obj.neurons(kNeuron).channel) ')'];
                text(.5*diff(xd)+xd(1), .05*(yd(2)-yd(1))+yd(1),n, 'Color', 'b', 'FontSize', 8)
                ylim(yd)
                %     text(.5*diff(xd)+xd(1), .05*(yd(2)-yd(1))+yd(1), num2str(kNeuron,2), 'Color', 'b', 'FontSize', 8)
            end
            
            set(gcf, 'currentaxes', ax(end))
            text(.5,.5, obj.neurons(kNeuron).exname)
            xlim([0 2])
            ylim([0 1])
            
            sciencestagram(gcf, 10, [10 10])
            
        end
        
        %% plot Rasters
        function h=plotRasters(obj, varargin)
            %             p.addOptional('aligningField', 'motionon')
            %             p.addOptional('sortField', [])
            %             p.addOptional('frozen', false)
            %             p.addOptional('plotMode', 'image')
            %             p.addOptional('window', [-.2 1.5])
            %             p.addOptional('plotPSTH', false)
            %             p.addOptional('indices', []);
            p=inputParser();
            p.addOptional('aligningField', 'motionon')
            p.addOptional('sortField', [])
            p.addOptional('frozen', false)
            p.addOptional('plotMode', 'image')
            p.addOptional('window', [-.2 1.5])
            p.addOptional('plotPSTH', false)
            p.addOptional('indices', []);
            p.parse(varargin{:});
            win=p.Results.window;
            h=gcf;
            nNeurons=numel(obj.neurons);
            sx=5;
            sy=5;
            ax=tight_subplot(sy,sx,.002, .002);
            for k=1:numel(ax)
                set(gcf, 'currentaxes', ax(k))
                axis off
            end
            
            for kNeuron=1:nNeurons
                set(gcf, 'currentaxes', ax(kNeuron))
                id=obj.neurons(kNeuron).trialIndex;
                if p.Results.frozen
                    [~, ii]=max(obj.stim.trialCnt);
                    iix=find(obj.stim.trialId==ii);
                    id=intersect(id, iix);
                end
                aligningField=[obj.neurons(kNeuron).stim.timing(id).(p.Results.aligningField)];
                
                if ~isempty(p.Results.sortField)
                    sortField=[obj.neurons(kNeuron).stim.timing(id).(p.Results.sortField)];
                    [~, ii]=sort(sortField-aligningField, 'ascend');
                else
                    ii=1:numel(aligningField);
                end
                ev=aligningField+[obj.neurons(kNeuron).stim.timing(id).plxstart];
                
                spcnt=binSpTimes(obj.neurons(kNeuron).spikeTimes, ev(ii), win, 1e-3);
                switch p.Results.plotMode
                    case 'image'
                        if ~isempty(p.Results.indices)
                            if iscell(p.Results.indices)
                                for k=1:numel(p.Results.indices)
                                    imagesc(p.Results.indices{k}, 1:size(spcnt,2), -spcnt(p.Results.indices{k},:)+1, [0 .1]); colormap gray
                                    hold on
                                end
                            else
                                imagesc(p.Results.indices, 1:size(spcnt,2), -spcnt(p.Results.indices,:)+1, [0 .1]); colormap gray
                            end
                        else
                            imagesc(-spcnt+1, [0 .1]); colormap gray
                        end
                    case 'point'
                        if ~isempty(p.Results.indices)
                            if iscell(p.Results.indices)
                                for k=1:numel(p.Results.indices)
                                    [i,j]=find(spcnt(p.Results.indices{k},:));
                                    plot(j,p.Results.indices{k}(i),'.k', 'MarkerSize', 1); hold on
                                end
                            else
                                [i,j]=find(spcnt(p.Results.indices,:));
                                plot(j,p.Results.indices(i),'.k', 'MarkerSize', 1); hold on
                            end
                        else
                            [i,j]=find(spcnt);
                            plot(j,i,'.k', 'MarkerSize', 1); hold on
                            
                        end
                    case 'line'
                end
                axis on
                hold on
                set(gca, 'XTick', '', 'YTick', '')
                axis tight
                yd=ylim;
                xd=xlim;
                n=[num2str(kNeuron,2) ' (' num2str( obj.neurons(kNeuron).channel) ')'];
                text(.5*diff(xd)+xd(1), .05*(yd(2)-yd(1))+yd(1),n, 'Color', 'b', 'FontSize', 8)
                ylim(yd)
                if p.Results.plotPSTH
                    ax2=axes('Position', get(gca, 'Position'), 'Color', 'none');
                    
                    plot(ax2,smooth(mean(spcnt),10), 'r')
                    xlim(xd)
                    yy=ylim;
                    yy(2)=yy(2)*2;
                    ylim(yy)
                    axis off
                end
                %     text(.5*diff(xd)+xd(1), .05*(yd(2)-yd(1))+yd(1), num2str(kNeuron,2), 'Color', 'b', 'FontSize', 8)
            end
            
            set(gcf, 'currentaxes', ax(end))
            text(.5,.5, obj.neurons(kNeuron).exname)
            xlim([0 2])
            ylim([0 1])
            
        end
        
        %% plot trial
        function h = plotTrial(obj, kTrial)
            clrMT=[70 144 205]/255;
            clrLIP=[85, 184, 74]/255;
            
            h=gca;
            hold all
            
            
            if nargin < 2
                kTrial = randi(numel(obj.trial));
            end
            
            if isfield(obj.trial(kTrial), 'eyepos')
                plot(1+zscore(obj.trial(kTrial).eyepos(:,1)), 'k');%(obj.trial(kTrial).eyepos(:,1)-mean(obj.trial(kTrial).eyepos(3e3:4e3,1)))/20, 'k'); hold on
                plot(1+zscore(obj.trial(kTrial).eyepos(:,2)), 'Color', .5*[1 1 1]);
                %                 plot(1+(obj.trial(kTrial).eyepos(:,2)-mean(obj.trial(kTrial).eyepos(3e3:4e3,2)))/20, 'Color', .5*[1 1 1]);
            end
            motion=zeros(obj.trial(kTrial).duration,1);
            pdur=150; %mean(diff(trial.pulseon));
            for kPulse=1:7
                motion(obj.trial(kTrial).pulseon(kPulse):obj.trial(kTrial).pulseon(kPulse)+pdur+16)=obj.trial(kTrial).pulses(kPulse);
            end
            contrast=zeros(obj.trial(kTrial).duration,1);
            targets=zeros(obj.trial(kTrial).duration,1);
            fixation=zeros(obj.trial(kTrial).duration,1);
            contrast(obj.trial(kTrial).motionon:obj.trial(kTrial).motionoff)=1;
            targets(obj.trial(kTrial).targson:obj.trial(kTrial).targsoff)=1;
            fixation(obj.trial(kTrial).fpon:obj.trial(kTrial).gosignal)=1;
            plot(2+.7*fixation, 'k');
            plot(3+.7*contrast, 'k')
            plot(4+.7*motion*3, 'k');
            plot(5+.7*targets, 'k')
            
            flabel{1}='eye xy';
            flabel{2}='fixation';
            flabel{3}='contrast';
            flabel{4}='direction';
            flabel{5}='targets';
            
            fctr=6;
            fields=setdiff(findFile(fieldnames(obj.trial(kTrial)), 'neuron'), findFile(fieldnames(obj.trial(kTrial)), 'FR_'));
            for kField=1:numel(fields)
                tmpx=obj.trial(kTrial).(fields{kField});
                tmpy=fctr*ones(size(tmpx));
                if any(strfind(fields{kField}, 'LIP'))
                    plot([tmpx(:) tmpx(:)]', [tmpy(:) tmpy(:)+1]', '-', 'Color', clrLIP)
                elseif any(strfind(fields{kField}, 'MT'))
                    plot([tmpx(:) tmpx(:)]', [tmpy(:) tmpy(:)+1]', '-', 'Color', clrMT)
                else
                    plot([tmpx(:) tmpx(:)]', [tmpy(:) tmpy(:)+1]', '-k')
                end
                flabel{fctr} = fields{kField};
                fctr = fctr + 1;
                
            end
            set(gca, 'Ytick', [1:5 (6:fctr-1)+.5], 'YtickLabel', flabel)
            
            
            
            title(sprintf('Trial %d', kTrial))
        end

        
        function plotChoPSTH(obj, n1, event, varargin)
            % compute a cross-covariogram for spikes aligned to an event
            % computeEventRelatedCovariogram(trial, neuron1, neuron2, eventName, varargin)
            
            p=inputParser();
            p.addOptional('lag', 200)
            p.addOptional('smoothing', 20)
            p.addOptional('window', [-100 1e3])
            p.parse(varargin{:})
            
            if isempty(obj.binSize)
                obj.binSize=10;
            end
            
            smwin=p.Results.smoothing/obj.binSize;
            
            n=neuroGLM('binSize', obj.binSize);
            
            neuronNames=findFile(fieldnames(obj.trial), 'neuron');
            assert(any(strcmp(neuronNames, n1)), 'Neuron names must be fields of trial struct')
            
            y1=full(n.getBinnedSpikeTrain(obj.trial, n1, 1:numel(obj.trial)));
            if smwin>0
                y1=smooth(y1, smwin);
            end
            
            ev=find(n.getBinnedSpikeTrain(obj.trial, event, 1:numel(obj.trial)));
            win=p.Results.window/obj.binSize;
            lag=max(abs(win));
            iter=0;
            maxIter=100;
            while any(diff(ev)<lag) && iter < maxIter
                ev(diff(ev)<lag)=[];
                iter=iter+1;
            end
%             fprintf('exited loop at iter %d. %d events left\n', iter, numel(ev))
            
            
            %             coh=arrayfun(@(x) sum(x.pulses), obj.trial);
            cho=arrayfun(@(x) x.choice, obj.trial);
            
            
            psthTime=(win(1):win(2))*obj.binSize;
            nPsthBins=numel(psthTime);
            rtrue=binContinuous(y1, ev, win);
            %             dpcell=dprime(nanmean(rtrue(:,psthTime>0 & psthTime < 1e3),2),cho);
            %             if dpcell<0 && choFlip
            %                 cho=~cho;
            %             end
            
            
            rbar=nan(nPsthBins,2);
            rbar(:,1)=nanmean(rtrue(cho,:));
            rbar(:,2)=nanmean(rtrue(~cho,:));
            
            
            set(gcf, 'DefaultAxesColorOrder', cbrewer('jake', 'rdbu', size(rbar,2)))
            plot(psthTime, rbar/obj.binSize*1e3)
            hold on
            plot(psthTime, nanmean(rtrue)/obj.binSize*1e3, 'k', 'Linewidth', 2)
            title(n1)
            xlabel('ms')
            ylabel('sp s^{-1}')
            
        end
        
        
        
        function plotCohPSTH(obj, n1, event, varargin)
            % compute a coherence-sorted PSTH
            % plotCohPSTH(obj, n1, event, varargin)
            
            p=inputParser();
            p.addOptional('lag', 200)
            p.addOptional('smoothing', 100)
            p.addOptional('window', [-100 1e3])
            p.parse(varargin{:})
            
            if isempty(obj.binSize)
                obj.binSize=10;
            end
            
            smwin=p.Results.smoothing/obj.binSize;
            
            n=neuroGLM('binSize', obj.binSize);
            
            neuronNames=findFile(fieldnames(obj.trial), 'neuron');
            assert(any(strcmp(neuronNames, n1)), 'Neuron names must be fields of trial struct')
            
            y1=full(n.getBinnedSpikeTrain(obj.trial, n1, 1:numel(obj.trial)));
            if smwin>0
                y1=smooth(y1, smwin);
            end
            
            ev=find(n.getBinnedSpikeTrain(obj.trial, event, 1:numel(obj.trial)));
            win=p.Results.window/obj.binSize;
            lag=max(abs(win));
            iter=0;
            maxIter=100;
            while any(diff(ev)<lag) && iter < maxIter
                ev(diff(ev)<lag)=[];
                iter=iter+1;
            end
%             fprintf('exited loop at iter %d. %d events left\n', iter, numel(ev))
            
            
            coh=arrayfun(@(x) sum(x.pulses), obj.trial);
            cho=arrayfun(@(x) x.choice, obj.trial);
            
            
            psthTime=(win(1):win(2))*obj.binSize;
            nPsthBins=numel(psthTime);
            rtrue=binContinuous(y1, ev, win);
            %             dpcell=dprime(nanmean(rtrue(:,psthTime>0 & psthTime < 1e3),2),cho);
            %             if dpcell<0 && choFlip
            %                 cho=~cho;
            %             end
            
            abscoh=abs(coh);
            abscoh(abscoh==0)=.01;
            scoh=sign(cho-.5).*abscoh; % aligned to choice now
            q=quantile(abs(coh), [.3 .6 .9 1]);
            binEdges = sort([-q 0 q]);
            nCohs=numel(binEdges)-1;
            binId=cell2mat(arrayfun(@(x,y) scoh>=x & scoh<=y, binEdges(1:end-1), binEdges(2:end), 'UniformOutput', false));
            
            
            rbar=nan(nPsthBins, nCohs);
            for k=1:nCohs
                rbar(:,k)=nanmean(rtrue(binId(:,k),:));
            end
            

            set(gcf, 'DefaultAxesColorOrder', cbrewer('jake', 'rdbu', size(rbar,2)))
            plot(psthTime, rbar/obj.binSize*1e3)
            hold on
            plot(psthTime, nanmean(rtrue)/obj.binSize*1e3, 'k', 'Linewidth', 2)
            title(n1)
            xlabel('ms')
            
        end % plot Coh PSTH
        
    end % methods
    
end

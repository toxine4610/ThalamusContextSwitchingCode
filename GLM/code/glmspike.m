% glmspike create a generalized linear model object for spike data
% This is going to get annoyingly large
classdef glmspike < neuroGLM
    
    properties
        description
        model
        modelfit
        designOptions
        fitNeuron
    end
    
    methods
        %% Constructor
        function obj = glmspike(desc)
            if nargin==1
                obj.description=desc;
            end
            % setup defaults
            obj.designOptions.normalizeBases     = true;
            obj.designOptions.orthogonalizeBases =false;
            obj.designOptions.unitOfTime         ='ms';
            obj.designOptions.binSize            = 10;
            obj.designOptions.includeSpikes 	= false;
            obj.designOptions.includeMT     	= false;
            obj.designOptions.includeLIP    	= false;
            obj.designOptions.MTbasis       	= 0;
            obj.designOptions.LIPbasis      	= 0;
            obj.designOptions.includeHistory	= false;
            obj.designOptions.historyField  	= [];
            obj.designOptions.boxcarPulse   	= false;
            obj.designOptions.boxcarContrast    = true;
            obj.designOptions.contrastPulse 	= false;
            obj.designOptions.includeContrast	= true;
            obj.designOptions.includeChoice     = false;
            obj.designOptions.separatePulses    = false;
            obj.designOptions.pulseBasis        = [];
            obj.designOptions.includeMotionOnset=false;
            obj.designOptions.includeMotionOffset= false;
            obj.designOptions.includeTargets    = false;
            obj.designOptions.includeFpoff      = false;
            obj.designOptions.includeFpon       = false;
            obj.designOptions.includeFpentered  = false;
            obj.designOptions.instantaneousCoupling = false;
            obj.designOptions.pulseType         = 'MT';
            obj.designOptions.pulseMultiplier   = [];
            obj.designOptions.CouplingEpoch     = [];
            obj.designOptions.MTCouplingWindow  = 800;
            obj.designOptions.excludeMotion=false;
            obj.designOptions.MTsim=false;
            obj.designOptions.trialDC = false;
            obj.designOptions.choiceDependentCoupling=false;
            obj.designOptions.couplingWindow =false;
            
            obj.model.regressionMode = 'RIDGEFIXED';
            
        end
        
        %% Setup design options
        function setDesignOptions(obj, varargin)
            p = inputParser();
            p.addOptional('normalizeBases',     obj.designOptions.normalizeBases);
            p.addOptional('orthogonalizeBases', obj.designOptions.orthogonalizeBases);
            p.addOptional('unitOfTime',         obj.designOptions.unitOfTime);
            p.addOptional('binSize',            obj.designOptions.binSize);
            p.addOptional('includeSpikes',      obj.designOptions.includeSpikes);
            p.addOptional('includeMT',          obj.designOptions.includeMT);
            p.addOptional('includeLIP',         obj.designOptions.includeLIP);
            p.addOptional('MTbasis',            obj.designOptions.MTbasis);
            p.addOptional('LIPbasis',           obj.designOptions.LIPbasis);
            p.addOptional('includeHistory',     obj.designOptions.includeHistory);
            p.addOptional('historyField',       obj.designOptions.historyField);
            p.addOptional('boxcarPulse',        obj.designOptions.boxcarPulse);
            p.addOptional('contrastPulse',      obj.designOptions.contrastPulse);
            p.addOptional('includeContrast',    obj.designOptions.includeContrast);
            p.addOptional('includeChoice',      obj.designOptions.includeChoice);
            p.addOptional('separatePulses',     obj.designOptions.separatePulses);
            p.addOptional('pulseBasis',         obj.designOptions.pulseBasis);
            p.addOptional('includeMotionOnset', obj.designOptions.includeMotionOnset);
            p.addOptional('includeMotionOffset',obj.designOptions.includeMotionOffset);
            p.addOptional('includeTargets',     obj.designOptions.includeTargets);
            p.addOptional('includeFpoff',       obj.designOptions.includeFpoff);
            p.addOptional('includeFpon',        obj.designOptions.includeFpon);
            p.addOptional('includeFpentered',   obj.designOptions.includeFpentered);
            p.addOptional('pulseType',          obj.designOptions.pulseType);
            p.addOptional('pulseMultiplier',    obj.designOptions.pulseMultiplier);
            p.addOptional('CouplingEpoch',      obj.designOptions.CouplingEpoch)
            p.addOptional('MTCouplingWindow',   obj.designOptions.MTCouplingWindow)
            p.addOptional('excludeMotion',      obj.designOptions.excludeMotion)
            p.addOptional('instantaneousCoupling', obj.designOptions.instantaneousCoupling)
            p.addOptional('MTsim',                  obj.designOptions.MTsim)
            p.addOptional('couplingWindow',         obj.designOptions.MTsim)
            p.addOptional('choiceDependentCoupling',obj.designOptions.choiceDependentCoupling)
            p.addOptional('trialDC',      obj.designOptions.trialDC)
            p.addOptional('boxcarContrast', obj.designOptions.boxcarContrast)
            p.parse(varargin{:});
            
            obj.designOptions=p.Results;
            
            assert(~isempty(obj.designOptions.binSize), 'binSize is empty')
            obj.binSize=obj.designOptions.binSize;
            obj.binfun = str2func(['@(t) (t==0) + ceil(t/' num2str(obj.binSize) ')']);
            
        end
        
        %% basis factory presets
        function bs = basisFactoryPresets(obj, basisName, duration, nBasis, nlScale)
            % Some Preset Basis styles
            % bs = basisFactoryPresets(basisName, duration, nBasis)
            % Possible basisName values:
            % history
            % spike
            % MTLIPcoupling
            % MTcoupling
            % LIPcoupling
            % MTmotion
            % LIPmotion
            % MTpulse
            % LIPpulse
            % Targets
            % Saccade
            if nargin ==1
                help glmspike/basisFactoryPresets
                if nargout>0
                    bs=[];
                end
                return
            end
            
            switch basisName
                case {'history', 'spike'}
                    if ~exist('duration', 'var'), duration=250;end
                    if ~exist('nBasis', 'var'), nBasis=8; end
                    if ~exist('nlScale', 'var'), nlScale=50; end
                    bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, nlScale);
                    % add delta functions to the basis for refractory period
                    nDelta=1;
                    bs.tr = [bs.tr(:,1:nDelta) bs.tr];
                    bs.B(1:nDelta,1:nDelta)=0; % first time bin is controlled by the delta functions
                    bs.B = [[eye(nDelta); zeros(size(bs.B,1)-nDelta,nDelta)] bs.B];
                    bs.edim = bs.edim + nDelta;
                    bs.nBases=bs.nBases+nDelta;
                    bs.centers=[1 bs.centers];
                case 'stim'
                    if ~exist('duration', 'var'), duration=[20 500];end
                    if ~exist('nBasis', 'var'), nBasis=8; end
                    if ~exist('nlScale', 'var'), nlScale=1e3; end
                    bs = basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, nlScale);
                    
                case {'MTLIPcoupling', 'InterAreaCoupling'}
                    if ~exist('duration', 'var'), duration=[20 800];end
                    if ~exist('nBasis', 'var'), nBasis=8; end
                    if ~exist('nlScale', 'var'), nlScale=250; end
                    bs = basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, nlScale);
                    
                case 'MTcoupling'
                    if ~exist('duration', 'var'), duration=[1 200];end
                    if ~exist('nBasis', 'var'), nBasis=8; end
                    bs = basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, nlScale);
                    
                case 'LIPcoupling'
                    if ~exist('duration', 'var'), duration=100;end
                    if ~exist('nBasis', 'var'), nBasis=8; end
                    bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, 2);
                    
                case 'LIPpulse'
                    if ~exist('duration', 'var'), duration=[50 1e3];end
                    if ~exist('nBasis', 'var'), nBasis=10; end
                    bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, 400);
                    
                case 'MTpulse'
                    if ~exist('duration', 'var'), duration=[50 400];end
                    if ~exist('nBasis', 'var'), nBasis=10; end
                    bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, 250);
                    
                case 'MTmotion'
                    if ~exist('duration', 'var'), duration=[0 1200];end
                    if ~exist('nBasis', 'var'), nBasis=20; end
                    bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, 500);
                    
                case 'LIPmotion'
                    if ~exist('duration', 'var'), duration=[50 1200];end
                    if ~exist('nBasis', 'var'), nBasis=12; end
                    bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, nBasis, obj.binfun, 500);
                case 'Targets'
                    if ~exist('duration', 'var'), duration=1500;end
                    if ~exist('nBasis', 'var'), nBasis=15; end
                    bs=basisFactory.makeSmoothTemporalBasis('raised cosine', duration, nBasis, obj.binfun);
                case {'Saccade', 'saccade', 'resp'}
                    if ~exist('duration', 'var'), duration=1000;end
                    if ~exist('nBasis', 'var'), nBasis=20; end
                    bs=basisFactory.makeSmoothTemporalBasis('raised cosine', duration, nBasis, obj.binfun);
                    
                otherwise
                    error('basis name not recognized')
                    
            end
            
            if obj.designOptions.normalizeBases
                bs.normalize;
            end
            
            if obj.designOptions.orthogonalizeBases
                bs.orthogonalize;
            end
        end
        
        %% buildDesign using boxcar for the different stimuli
        function buildDesignBoxcar(obj, trial, gammaStim, gammaHist, gammaCoupleSame, gammaCoupleInter)
            % build design specifications for neuron with different integration
            % parameters for the different types of regressors
            % buildDesignBoxcar(obj, trial, gammaStim, gammaHist, gammaCoupleSame, gammaCoupleInter, couplingId)
            % Inputs:
            %   trial      [mTrials x 1] trial struct array
            %   gammaStim        [1 x 1] integration window for stimulus
            %   gammaHist        [1 x 1] integration window for history filter
            %   gammaCoupleSame  [1 x 1] integration window for coupling filters
            %                            for same area neurons
            %   gammaCoupleInter [1 x 1] integration window for inter area
            %                            coupling filters
            % gammas entered in msec
            if nargin==1
                help glmspike/buildDesignBoxcar
                return
            end
            
            if ~exist('gammaStim', 'var')
                gammaStim = 500;
            end
            
            if ~exist('gammaHist', 'var')
                gammaHist = 250;
            end
            
            if ~exist('gammaCoupleSame', 'var')
                gammaCoupleSame = 250;
            end
            
            if ~exist('gammaCoupleInter', 'var')
                gammaCoupleInter = 500;
            end
            
            
            assert(~isempty(obj.fitNeuron), 'you must identify which neuron you are fitting')
            assert(any(strcmp(fieldnames(trial), obj.fitNeuron)), 'fitNeuron must exist in trial struct')
            
            % gather info about the neurons that could be covariates in
            % this dataset
            brainArea        = obj.fitNeuron(1:strfind(obj.fitNeuron, 'neuron')-1);
            neuronNames      = findFile(fieldnames(trial), 'neuron');
            sameAreaNeurons  = findFile(neuronNames, brainArea);
            otherAreaNeurons = setdiff(neuronNames, sameAreaNeurons);
            sameAreaNeurons  = setdiff(sameAreaNeurons, obj.fitNeuron);
            
            % --- regressor for pre-trial firing rate of the neuron
            if obj.designOptions.trialDC
                   stimHandle = @(trial) ones(obj.binfun(trial.duration), 1) * trial.(['FR_' obj.fitNeuron]);
                   B = Basis;
                   B.binfun = obj.binfun;
                   obj.addCovariate(trial, 'TrialDC', 'Trial DC', stimHandle, B);                
            end
            
            % -------------------------------------------------------------
            % --- Stimulus Covariates
            if gammaStim>0 % if any stimulus terms are included
                
                % --- Setup Basis functions for MT and LIP
                if strcmp(obj.fitNeuron(1:2), 'MT') % MT
                    stimBasis = basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', [30 gammaStim],10, obj.binfun, 100);
                    stimBasis.normalize;
                    
                    targBasis = basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', [30 gammaStim], 8, obj.binfun, 10);
                    targBasis.normalize;
                    
                else % LIP
                    stimBasis = obj.basisFactoryPresets('stim', [100 gammaStim], 10, 500); % longer integration time
                    targBasis = obj.basisFactoryPresets('Targets');
                end
                
                % add Targets
                obj.addCovariateTiming(trial, 'targson', 'targson', 'Targets Onset', targBasis)
                
                % --- Choice terms
                if obj.designOptions.includeChoice
                    
                    % if include choice is a number, use it as the offset
                    if abs(obj.designOptions.includeChoice)>1
                        offset=obj.designOptions.includeChoice;
                        nBases=ceil((abs(offset)+500)/70);
                    else
                        offset=-2500; %  ms before saccade onset
                        nBases=40;
                    end
                    
                    saccadeBasis = obj.basisFactoryPresets('saccade', abs(offset)+500, nBases);
                    obj.addCovariateTiming(trial, 'resp1', 'saccade', 'Saccade Choice 1', saccadeBasis, offset, @(trial) (trial.choice == 1));
                    obj.addCovariateTiming(trial, 'resp2', 'saccade', 'Saccade Choice 2', saccadeBasis, offset, @(trial) (trial.choice == 0));
                    
                else % just saccade kernels 500ms perisaccade
                    offset=-500; % 500 ms before saccade onset
                    saccadeBasis = obj.basisFactoryPresets('saccade', 1e3, 20);
                    obj.addCovariateTiming(trial, 'resp1', 'saccade', 'Saccade Choice 1', saccadeBasis, offset, @(trial) (trial.choice == 1));
                    obj.addCovariateTiming(trial, 'resp2', 'saccade', 'Saccade Choice 2', saccadeBasis, offset, @(trial) (trial.choice == 0));
                end
                
                
                %----------------------------------------------------------
                % Covariates related to the motion stimulus
                if ~obj.designOptions.excludeMotion % are directional terms included
                    
                    % --- use simulated MT as input in place of motion stimulus
                    if obj.designOptions.MTsim
                        
                        % --- include nondirectional contrast term as well
                        if obj.designOptions.includeContrast
                            if obj.designOptions.boxcarContrast
                                % --- add contrast term as a boxcar
                                obj.addCovariateBoxcar(trial, 'Motion', 'motionon', 'motionoff', [], 'Motion Contrast', stimBasis)
                            else
                                % --- add contrast using spatially averaged contrast
                                if obj.binSize==1
                                    stimHandle = @(trial) trial.contrast;
                                else
                                    stimHandle = @(trial) basisFactory.binDataVector(trial.contrast/obj.binSize, obj.binSize);
                                end
                                obj.addCovariate(trial, 'Motion', 'Motion Contrast', stimHandle, stimBasis);
                            end
                        end
                        
                        
                        % add MT pref
                        if obj.binSize==1
                            prefHandle = @(trial) trial.MTpref;
                        else
                            prefHandle = @(trial) basisFactory.binDataVector(trial.MTpref/obj.binSize, obj.binSize);
                        end
                        obj.addCovariate(trial,'MTpref', 'MT pref', prefHandle, stimBasis);
                        
                        % add MT anti
                        if obj.binSize==1
                            antiHandle = @(trial) trial.MTanti;
                        else
                            antiHandle = @(trial) basisFactory.binDataVector(trial.MTanti/obj.binSize, obj.binSize);
                        end
                        obj.addCovariate(trial,'MTanti', 'MT anti', antiHandle, stimBasis);
                        
                        % --- or, construct pulse input as boxcars
                    elseif obj.designOptions.boxcarPulse
                        
                        if obj.designOptions.boxcarContrast
                            % --- add contrast term as a boxcar
                            obj.addCovariateBoxcar(trial, 'Motion', 'motionon', 'motionoff', [], 'Motion Contrast', stimBasis)
                        else
                            % --- add contrast using spatially averaged contrast
                            if obj.binSize==1
                                stimHandle = @(trial) trial.contrast;
                            else
                                stimHandle = @(trial) basisFactory.binDataVector(trial.contrast/obj.binSize, obj.binSize);
                            end
                            obj.addCovariate(trial, 'Motion', 'Motion Contrast', stimHandle, stimBasis);
                        end
                        
                        % --- add motion pulses
                        obj.addCovariateBoxcar(trial, 'Pulse', 'pulseon', 'pulseoff', 'pulses', 'Motion Pulses', stimBasis);
                        
                        % --- otherwise, pulses are impulses
                    else
                        motionOnsetBasis=obj.basisFactoryPresets('stim', [50 1100], 30, 500);
                        obj.addCovariateTiming(trial, 'motionon', 'motionon', 'Motion Stimulus Onset', motionOnsetBasis);
                        obj.addCovariateTiming(trial, 'motionoff', 'motionoff', 'Motion Stimulus Offset', stimBasis);
                        stimHandle = @(trial) sparse(obj.binfun(trial.pulseon), 1, trial.pulses, obj.binfun(trial.duration), 1);
                        obj.addCovariate(trial,'Pulse', 'Motion Pulses', stimHandle,stimBasis);
                    end % Selection of stimulus input (if MTsim)
                    
                    
                    
                    if obj.designOptions.includeMotionOnset
                        bs=obj.basisFactoryPresets('stim', [50 1100], 30, 500);
                        obj.addCovariateTiming(trial, 'Onset', 'motionon', 'Motion Onset', bs, 0);
                    end
                    
                    if obj.designOptions.includeMotionOffset
                        bs=obj.basisFactoryPresets('stim', [50 500], 10, 500);
                        obj.addCovariateTiming(trial, 'Onset', 'motionon', 'Motion Onset', bs, 0);
                    end
                    
                else % exclude directional terms
                    % --- include nondirectional contrast term as well
                    if obj.designOptions.includeContrast
                        if obj.designOptions.boxcarContrast
                            % --- add contrast term as a boxcar
                            obj.addCovariateBoxcar(trial, 'Motion', 'motionon', 'motionoff', [], 'Motion Contrast', stimBasis)
                        else
                            % --- add contrast using spatially averaged contrast
                            if obj.binSize==1
                                stimHandle = @(trial) trial.contrast;
                            else
                                stimHandle = @(trial) basisFactory.binDataVector(trial.contrast/obj.binSize, obj.binSize);
                            end
                            obj.addCovariate(trial, 'Motion', 'Motion Contrast', stimHandle, stimBasis);
                        end
                    end
                    
                end % if directional terms are included
                
            end % if any stimulus terms are included
            
            % -------------------------------------------------------------
            % --- History Filter
            if gammaHist>0
                historyBasis=obj.basisFactoryPresets('history', gammaHist);
                obj.addCovariateSpiketrain(trial, 'History', obj.fitNeuron, obj.fitNeuron, historyBasis);
            end
            
            % -------------------------------------------------------------
            % --- Inter-Area Coupling
            
            nOtherNeurons=numel(otherAreaNeurons);
            
            couplingId=1:nOtherNeurons;
            
            
            if nOtherNeurons>0 && gammaCoupleInter>0
                couplingBasis=obj.basisFactoryPresets('InterAreaCoupling', [0 gammaCoupleInter], 8, 100);
                
                for kNeuron=couplingId(:)'
                    
                    % --- are spikes counts or times
                    if numel(trial(1).(otherAreaNeurons{kNeuron}))==trial(1).duration % counts?

                        if obj.binSize==1
                            stimHandle = @(trial) trial.(otherAreaNeurons{kNeuron});
                        else
                            stimHandle = @(trial) basisFactory.binDataVector(mean(trial.(otherAreaNeurons{kNeuron}), 2)/obj.binSize, obj.binSize);
                        end
                        
                        % --- check if coupling depends on the choice
                        if obj.designOptions.choiceDependentCoupling
                            obj.addCovariateTiming(trial, sprintf('InterAreaCoupling%dChoice1', kNeuron), otherAreaNeurons{kNeuron}, [otherAreaNeurons{kNeuron} 'choice1'], couplingBasis, 1, @(trial) (trial.choice == 1));
                            obj.addCovariateTiming(trial, sprintf('InterAreaCoupling%dChoice2', kNeuron), otherAreaNeurons{kNeuron}, [otherAreaNeurons{kNeuron} 'choice2'], couplingBasis, 1, @(trial) (trial.choice == 0));
                        else
                            obj.addCovariate(trial, sprintf('InterAreaCoupling%d', kNeuron), otherAreaNeurons{kNeuron}, stimHandle, couplingBasis);
                        end
                        
                    else % spikes are saved as times
                        
                        if obj.designOptions.instantaneousCoupling % no offset in simultaneous spikes (instantaneous coupling)
                            
                            if obj.designOptions.choiceDependentCoupling
                                obj.addCovariateTiming(trial, sprintf('InterAreaCoupling%dChoice1', kNeuron), otherAreaNeurons{kNeuron}, [otherAreaNeurons{kNeuron} 'choice1'], couplingBasis, 0, @(trial) (trial.choice == 1));
                                obj.addCovariateTiming(trial, sprintf('InterAreaCoupling%dChoice2', kNeuron), otherAreaNeurons{kNeuron}, [otherAreaNeurons{kNeuron} 'choice2'], couplingBasis, 0, @(trial) (trial.choice == 0));
                            else
                                obj.addCovariateTiming(trial, sprintf('InterAreaCoupling%d', kNeuron), otherAreaNeurons{kNeuron}, otherAreaNeurons{kNeuron}, couplingBasis);
                            end
                            
                        else % only causal interactions (coupling lagged 1 bin)
                            if obj.designOptions.choiceDependentCoupling
                                obj.addCovariateTiming(trial, sprintf('InterAreaCoupling%dChoice1', kNeuron), otherAreaNeurons{kNeuron}, [otherAreaNeurons{kNeuron} 'choice1'], couplingBasis, 1, @(trial) (trial.choice == 1));
                                obj.addCovariateTiming(trial, sprintf('InterAreaCoupling%dChoice2', kNeuron), otherAreaNeurons{kNeuron}, [otherAreaNeurons{kNeuron} 'choice2'], couplingBasis, 1, @(trial) (trial.choice == 0));
                            else
                                obj.addCovariateSpiketrain(trial, sprintf('InterAreaCoupling%d', kNeuron), otherAreaNeurons{kNeuron}, otherAreaNeurons{kNeuron}, couplingBasis);
                            end
                        end % instantaneous coupling
                        
                    end % if spikes are times or counts
                    
                end % loop over coupled neurons
                
            end % if inter-areal coupling terms should exist
            
            
            % -------------------------------------------------------------
            % --- IntraAreaCoupling
            
            nSameNeurons=numel(sameAreaNeurons);
            
            couplingId=1:nSameNeurons;
            
            
            if nSameNeurons>0 && gammaCoupleSame>0
                couplingBasis=obj.basisFactoryPresets('InterAreaCoupling', [0 gammaCoupleSame], 8, 100);
                
                for kNeuron=couplingId(:)'
                    
                    % --- are spikes counts or times
                    if numel(trial(1).(sameAreaNeurons{kNeuron}))==trial(1).duration

                        if obj.binSize==1
                            stimHandle = @(trial) trial.(sameAreaNeurons{kNeuron});
                        else
                            stimHandle = @(trial) basisFactory.binDataVector(mean(trial.(otherAreaNeurons{kNeuron}), 2)/obj.binSize, obj.binSize);
                        end
                        
                        if obj.designOptions.choiceDependentCoupling
                            obj.addCovariateTiming(trial, sprintf('SameAreaCoupling%dChoice1', kNeuron), sameAreaNeurons{kNeuron}, [sameAreaNeurons{kNeuron} 'choice1'], couplingBasis, 1, @(trial) (trial.choice == 1));
                            obj.addCovariateTiming(trial, sprintf('SameAreaCoupling%dChoice2', kNeuron), sameAreaNeurons{kNeuron}, [sameAreaNeurons{kNeuron} 'choice2'], couplingBasis, 1, @(trial) (trial.choice == 0));
                        else
                            obj.addCovariate(trial, sprintf('SameAreaCoupling%d', kNeuron), sameAreaNeurons{kNeuron}, stimHandle, couplingBasis);
                        end
                        
                    else
                        if obj.designOptions.instantaneousCoupling
                            
                            if obj.designOptions.choiceDependentCoupling
                                obj.addCovariateTiming(trial, sprintf('SameAreaCoupling%dChoice1', kNeuron), sameAreaNeurons{kNeuron}, [sameAreaNeurons{kNeuron} 'choice1'], couplingBasis, 0, @(trial) (trial.choice == 1));
                                obj.addCovariateTiming(trial, sprintf('SameAreaCoupling%dChoice2', kNeuron), sameAreaNeurons{kNeuron}, [sameAreaNeurons{kNeuron} 'choice2'], couplingBasis, 0, @(trial) (trial.choice == 0));
                            else
                                obj.addCovariateTiming(trial, sprintf('SameAreaCoupling%d', kNeuron), sameAreaNeurons{kNeuron}, sameAreaNeurons{kNeuron}, couplingBasis);
                            end
                            
                        else
                            if obj.designOptions.choiceDependentCoupling
                                obj.addCovariateTiming(trial, sprintf('SameAreaCoupling%dChoice1', kNeuron), sameAreaNeurons{kNeuron}, [sameAreaNeurons{kNeuron} 'choice1'], couplingBasis, 1, @(trial) (trial.choice == 1));
                                obj.addCovariateTiming(trial, sprintf('SameAreaCoupling%dChoice2', kNeuron), sameAreaNeurons{kNeuron}, [sameAreaNeurons{kNeuron} 'choice2'], couplingBasis, 1, @(trial) (trial.choice == 0));
                            else
                                obj.addCovariateSpiketrain(trial, sprintf('SameAreaCoupling%d', kNeuron), sameAreaNeurons{kNeuron}, sameAreaNeurons{kNeuron}, couplingBasis);
                            end
                        end % if instantaneous coupling
                        
                    end % if counts or spikes
                    
                end % loop over coupled neurons
                
            end % if include within area coupling
            
            
        end % buildDesignBoxcar
        
        %% fit the model with a fixed ridge parameter and calculate AIC
        function M=fitRidgeFixed(obj, trial, rho)
            % fit poisson regression on full dataset with fixed ridge
            % M=fitRidgeFixed(obj, trial, rho)
            % parameter and calculate AIC
            
            if nargin==1
                help glmspike/fitRidgeFixed
                if nargout>0
                    M=[];
                end
                return
            end
            
            assert(~isempty(obj.dm.X), 'compile design matrix first!')
            
            if ~exist('rho', 'var')
                rho=.2;
            end
            
            obj.model.regressionMode='RIDGEFIXED';
            
            y=obj.getBinnedSpikeTrain(trial, obj.fitNeuron, obj.dm.trialIndices);
            
            M=regression.doRegressionPoisson(obj.dm.X, y, obj, [], obj.binSize/1e3, rho);
            
            lambda=M.fnlin(obj.dm.X*M.khat)*M.dt;
            M.logli=logliPoisson(lambda, y);
            M.df=obj.edim;
            M.AIC=2*obj.edim - 2*M.logli;
            
        end
        
        %% fit the model with cross validation
        function fitCV(obj,trial, nFolds, rho, colInds)
            % fitCV(obj,trial, nFolds)
            if nargin==1
                help fitCV
                return
            end
            
            if ~exist('nFolds', 'var')
                nFolds=5;
            end
            
            if ~exist('rho', 'var')
                rho=.2;
            end
            
            if ~exist('colInds', 'var')
                colInds=1:obj.edim;
            end
            
            y=obj.getBinnedSpikeTrain(trial, obj.fitNeuron, obj.dm.trialIndices);
            
            % cross validate it
            nFitTrials=numel(obj.dm.trialIndices);
            
            % get trial indices
            trialStartT  = [1 cumsum(obj.binfun([trial(obj.dm.trialIndices).duration]))];
            trialSampleIndices = cell(nFitTrials,1);
            
            obj.model.trialSampleIndices=trialSampleIndices;
            
            for k = 1:numel(obj.dm.trialIndices)
                trialSampleIndices{k} = trialStartT(k):trialStartT(k+1);
            end
            
            if nFolds == 1
                xvFolds = regression.xvalidationIdx(nFitTrials, 10, true, false);
            else
                xvFolds = regression.xvalidationIdx(nFitTrials, nFolds, false, true);
            end
            
            
            for kFold = 1:nFolds
                ndx = [trialSampleIndices{xvFolds{kFold,1}}];
                obj.modelfit(kFold).regressionMode=obj.model.regressionMode;
                obj.modelfit(kFold).trainingIndices=obj.dm.trialIndices(xvFolds{kFold,1});
                obj.modelfit(kFold).testIndices=obj.dm.trialIndices(xvFolds{kFold,2});
                obj.modelfit(kFold).dt=obj.binSize/1e3;
                
                tmp=regression.doRegressionPoisson(obj.dm.X, y, obj, ndx, obj.modelfit(kFold).dt, rho, colInds);
                obj.modelfit(kFold).khat=tmp.khat;
                obj.modelfit(kFold).fnlin=tmp.fnlin;
                obj.modelfit(kFold).SDebars=tmp.SDebars;
                if isfield(tmp, 'rho')
                    obj.modelfit(kFold).rho=tmp.rho;
                end
            end
        end
        
        %% predict spike rate
        function [rhat, trialStarts, trialEnds] = predictSpikeRate(obj, M, trial, varargin)
            % predict spike rate fiven stimulus covariates
            % [rhat, trialStarts, trialEnds] = predictSpikeRate(M, trial, varargin)
            % Optional Inputs
            %   'trialIndices',     default = dm.trialIndices
            %   'history',          default = []
            %   'includeHistory',   default = false
            %   'nlinfun',          default = [] (must be a function handle)
            %   'forceChoice',      default = 0 (1 or 2 forces all choices to be coded
            %                                 as specified)
            if nargin==1
                help predictSpikeRate
                return
            end
            p = inputParser();
            
            if ~isfield(obj.dm, 'trialIndices') || isempty(obj.dm.trialIndices)
                trialInds=1:numel(trial);
            else
                trialInds=obj.dm.trialIndices;
            end
            p.addOptional('trialIndices', trialInds);
            p.addOptional('includeHistory', true);
            p.addOptional('forceChoice', 0);
            p.addOptional('forceDirectionPreference', 0);
            p.addOptional('spiking', false);
            p.parse(varargin{:});
            
            
            ws=obj.combineWeights(M.khat);
            [X, trialStarts, trialEnds] = obj.getCovariateTiming(trial, p.Results.trialIndices);
            Covariates = {obj.covar.label};
            weightLabels=fieldnames(ws);
            pulseKernels=find(strcmp(weightLabels, 'Pulse'));
            
            if p.Results.forceDirectionPreference==1
                for kCov=1:numel(pulseKernels)
                    if sum(ws.(weightLabels{pulseKernels(kCov)}).data)<0
                        ws.(weightLabels{pulseKernels(kCov)}).data=-ws.(weightLabels{pulseKernels(kCov)}).data;
                    end
                end
            elseif p.Results.forceDirectionPreference==2
                for kCov=1:numel(pulseKernels)
                    if sum(ws.(weightLabels{pulseKernels(kCov)}).data)>0
                        ws.(weightLabels{pulseKernels(kCov)}).data=-ws.(weightLabels{pulseKernels(kCov)}).data;
                    end
                end
                
            end
            
            
            history=find(strcmp({obj.covar.desc}, obj.fitNeuron) | strcmp({obj.covar.desc}, 'Spike Hisory') | strcmp({obj.covar.desc}, 'Spike History'),1);
            nT = size(X,1);
            
            % This part just forces all choices to be 1 or 2
            if p.Results.forceChoice~=0
                switch p.Results.forceChoice
                    case 1
                        % motionChoice
                        covOrig = findFile(Covariates, {'motionChoice', '2'});
                        covReplace = findFile(Covariates, {'motionChoice', '1'});
                        if ~isempty(covOrig)
                            X(find(X(:, obj.idxmap.(covOrig{1}))), obj.idxmap.(covReplace{1})) = 1; %#ok<FNDSB>
                            X(:, obj.idxmap.(covOrig{1})) = 0;
                        end
                        % resp
                        covOrig     = findFile(Covariates, {'resp', '2'});
                        covReplace  = findFile(Covariates, {'resp', '1'});
                        if ~isempty(covOrig)
                            X(find(X(:, obj.idxmap.(covOrig{1}))), obj.idxmap.(covReplace{1})) = 1; %#ok<FNDSB>
                            X(:, obj.idxmap.(covOrig{1})) = 0;
                        end
                        
                    case 2
                        % motionChoice
                        covOrig     = findFile(Covariates, {'motionChoice', '1'});
                        covReplace  = findFile(Covariates, {'motionChoice', '2'});
                        if ~isempty(covOrig)
                            X(find(X(:, obj.idxmap.(covOrig{1}))), obj.idxmap.(covReplace{1})) = 1; %#ok<FNDSB>
                            X(:, obj.idxmap.(covOrig{1})) = 0;
                        end
                        % resp
                        covOrig     = findFile(Covariates, {'resp', '1'});
                        covReplace  = findFile(Covariates, {'resp', '2'});
                        if ~isempty(covOrig)
                            X(find(X(:, obj.idxmap.(covOrig{1}))), obj.idxmap.(covReplace{1})) = 1; %#ok<FNDSB>
                            X(:, obj.idxmap.(covOrig{1})) = 0;
                        end
                        
                    otherwise
                        disp('Choice must be 1 or 2')
                end
            end
            
            nTrials=numel(trialStarts);
            rhat = zeros(nT,1);
            if isfield(ws, 'bias')
                rhat = rhat + ws.bias;
                %                 rhat = ws.bias*ones(nT,1);
            end
            
            for kCov = 1:numel(Covariates)
                if ~isempty(history) && (strcmp(Covariates{kCov}, Covariates{history}) && ~p.Results.includeHistory)
                    disp('skipping history')
                    continue
                end
                
                if strcmp(Covariates{kCov}, 'bias')
                    continue
                end
                
                ii = obj.idxmap.(Covariates{kCov});
                
                for kTrial=1:nTrials
                    iix=trialStarts(kTrial):trialEnds(kTrial);
                    tmp = filter(ws.(Covariates{kCov}).data, 1, X(iix,ii));
                    rhat(iix) = rhat(iix) + tmp;
                end
            end
            
            % Still no good way to simulate with history
            % % %             rhat=rhat;
            % % %             rhat=(g.dm.X*model.wts);
            % %             lambda=rhat;
            % %             whist=flipud(ws.(Covariates{history}).data(1:10));
            % %
            % %             whist(end)=-1e3;
            % %             nh=numel(whist);
            % %             gy=gsmooth(y, 5);
            % %             for k=nh+1:numel(rhat)
            % %                 if k > 2e3
            % %                     rhat(k)=(model.fnlin(lambda(k) + rhat(k-nh:(k-1))'*whist)*(g.binSize/1e3))>0;
            % %                 else
            % %                     rhat(k)=(model.fnlin(lambda(k) )*(g.binSize/1e3));
            % %                 end
            % %
            % %                 if k>1e3
            % %                     ix=k-1000:k-1;
            % %                     plot(ix,gy(ix),ix,rhat(ix));
            % %                 end
            % %                 drawnow
            % %             end
            
            if ~isfield(M, 'dt')
                dt=1;
            else
                dt=M.dt;
            end
            rhat = M.fnlin(rhat)*dt; % rate
            if any(isnan(rhat))
                fprintf('rhat has nans\n')
            end
            if p.Results.spiking
                rhat=poissrnd(rhat);
            end
            
        end
        %% get weight struct
        function [wts,err]=getWeights(obj, M)
            % get a particular model weights
            %         [wts,err]=getWeights(obj, M)
            if nargin==1
                help glmspike/getWeights
                return
            end
            if nargout==1
            wts=obj.combineWeights(M.khat);
            end
            if nargout==2
                wts=obj.combineWeights(M.khat);
                err=obj.combineWeights(M.SDebars);
            end
            
            
            
        end
        
        %% plot weights
        function plotWeights(obj, M, plotSD)
%             plotWeights(obj, M, plotSD)

            if nargin <2
                help glmspike/plotWeights
                return
            end
            
            if ~exist('plotSD', 'var')
                plotSD=false;
            end
            if ~isstruct(M)
                wts=obj.combineWeights(M);
                plotSD=false;
            else
                [wts,sd]=obj.getWeights(M);
            end
            f=fieldnames(wts);
            f(strcmp(f, 'bias'))=[];
            tstr={obj.covar.desc};
            n=numel(f);
            sx=ceil(sqrt(n));
            figure();
            for k=1:n
                subplot(sx,sx,k)
                if plotSD
                    errorbarFill(wts.(f{k}).tr, wts.(f{k}).data, sd.(f{k}).data, 'k', 'FaceColor', .5*[1 1 1]); hold on
                end
                plot(wts.(f{k}).tr, wts.(f{k}).data, 'k'); hold on
                axis tight
                plot(xlim, [0 0], 'k')
                title(tstr{k})
            end
            
        end
        
        %% save
        function fs=save(obj, modelDir, Tag)
            % save glmspike object
            % save(obj, modelDir, Tag)
            if nargin==1
                help glmspike/save
                return
            end
            
            if ~exist('Tag', 'var')
                Tag='';
            end
            % remove design matrix
            obj.dm.X=[];
            disp('saving')
            fname=fullfile(modelDir, [Tag, obj.fitNeuron, obj.description]);
            try
                save(fname, '-v7', 'obj')
            catch
                save(fname, '-v7.3', 'obj')
            end
            disp('done')
            if nargout>0
                fs=fname;
            end
        end
        
        %% plot Coupling
        function ff=plotCoupling(obj, M)
            if nargin < 2
                M=obj.modelfit;
            end
            covariateDescriptions={obj.covar(:).desc};
            [couplingFields, covarIds]=findFile(covariateDescriptions, 'neuron');
            couplingLabels={obj.covar(covarIds).label};
            for kFold=1:numel(M)
                ws(kFold)=obj.getWeights(M(kFold));
            end
            nCov=numel(couplingLabels);
            if isempty(couplingLabels)
                return
            end
            
            cf=cell(nCov,1);
            for kCov=1:nCov
                subplot(nCov,1,kCov)
                data=cell2mat(arrayfun(@(x) x.data, [ws.(couplingLabels{kCov})], 'uniformoutput', false));
                m=mean(data,2);
                if strcmp((couplingLabels{kCov}), 'History')
                    cf{kCov}=nan;
                else
                    cf{kCov}=m;
                end
%                 if sign(sum(m))<0
%                     m=-m;
%                 end
                errorbar(ws(1).(couplingLabels{kCov}).tr, m, std(data,[],2)); hold on
                title(couplingFields{kCov})
                
            end
            
            if nargout==1
                ff=cf;
            end
            
            
        end
        
        %% plot motion and saccade
        function varargout = plotMotion(obj, M)
            if nargin < 2
                M=obj.modelfit;
            end
            covariateDescriptions={obj.covar(:).desc};
            [motionFields, covarIds]=findFile(covariateDescriptions, 'Motion Choice');
            if isempty(motionFields)
                [motionFields, covarIds]=findFile(covariateDescriptions, 'Motion Onset');
            end
            if isempty(motionFields)
                [motionFields, covarIds]=findFile(covariateDescriptions, 'Motion Contrast');
            end
            motionLabels={obj.covar(covarIds).label};
            for kFold=1:numel(M)
                ws(kFold)=obj.getWeights(M(kFold));
            end
            nCov=numel(motionLabels);
            for kCov=1:nCov
                data=cell2mat(arrayfun(@(x) x.data, [ws.(motionLabels{kCov})], 'uniformoutput', false));
                m=mean(data,2);
                if sign(sum(m))<0
                    m=-m;
                end
                errorbar(ws(1).(motionLabels{kCov}).tr, m, std(data,[],2)); hold on
                title(motionFields{kCov})
            end
            
            if nargout > 0
               varargout{1} = m; 
            end
        end
        
        %% plot saccade response kernel
        function plotResp(obj, M)
            if nargin < 2
                M=obj.modelfit;
            end
            covariateDescriptions={obj.covar(:).desc};
            [saccadeFields, covarIds]=findFile(covariateDescriptions, 'Saccade Choice');
            if isempty(saccadeFields)
                [saccadeFields, covarIds]=findFile(covariateDescriptions, 'Saccade');
            end
            saccadeLabels={obj.covar(covarIds).label};
            for kFold=1:numel(M)
                ws(kFold)=obj.getWeights(M(kFold));
            end
            nCov=numel(saccadeLabels);
            for kCov=1:nCov
                data=cell2mat(arrayfun(@(x) x.data, [ws.(saccadeLabels{kCov})], 'uniformoutput', false));
                m=mean(data,2);
                if sign(sum(m))<0
                    m=-m;
                end
                errorbar(ws(1).(saccadeLabels{kCov}).tr, m, std(data,[],2)); hold on
                title(saccadeFields{kCov})
            end
            
            
        end
        
        %% plot motion pulse kernels
        function [pulseWeights, pulseErr]=plotPulse(obj, M)
            if nargin < 2
                M=obj.modelfit;
            end
            covariateDescriptions={obj.covar(:).desc};
            [pulseFields, covarIds]=findFile(covariateDescriptions, 'Pulses');
            
            if isempty(pulseFields)
                [pulseFields, covarIds]=findFile(covariateDescriptions, 'Motion Pulses');  
            end
            
            pulseLabels={obj.covar(covarIds).label};
            for kFold=1:numel(M)
                [ws(kFold), err(kFold)]=obj.getWeights(M(kFold));
            end
            nCov=numel(pulseLabels);
            for kCov=1:nCov
                data=cell2mat(arrayfun(@(x) x.data, [ws.(pulseLabels{kCov})], 'uniformoutput', false));
                errd=cell2mat(arrayfun(@(x) x.data, [err.(pulseLabels{kCov})], 'uniformoutput', false));
                m=mean(data,2);
                if sign(sum(m))<0
                    m=-m;
                end
                plot(ws(1).(pulseLabels{kCov}).tr, m); hold on
                s=mean(errd,2);
                errorbar(ws(1).(pulseLabels{kCov}).tr, m, s); hold on
                title(pulseFields{kCov})
            end
            
            if nargout >0
                pulseWeights=m;
                if nargout>1
                    pulseErr=s;
                end
            end
            
            
        end
        
        %% cross-validate predict rate
        function [lambda,trialRates]=cvPredictRate(obj, trial)
%            [lambda,trialRates]=cvPredictRate(obj, trial)
            trialRates=cell(numel(trial),1);
            
            for kFold=1:numel(obj.modelfit)
                [rhat, trialStarts, trialEnds]=obj.predictSpikeRate(obj.modelfit(kFold), trial, 'trialIndices', obj.modelfit(kFold).testIndices);
                for kTrial=1:numel(obj.modelfit(kFold).testIndices)
                    trialRates{obj.modelfit(kFold).testIndices(kTrial)}=rhat(trialStarts(kTrial):trialEnds(kTrial));
                    if isempty(trialRates{obj.modelfit(kFold).testIndices(kTrial)})
                        trialRates{obj.modelfit(kFold).testIndices(kTrial)}=nan(trial(obj.modelfit(kFold).testIndices(kTrial)).duration,1);
                    end
                end
            end
            
            emptyInds=find(cellfun(@isempty,trialRates));
            if ~isempty(emptyInds)
                [rhat, trialStarts, trialEnds]=obj.predictSpikeRate(obj.modelfit(kFold), trial, 'trialIndices', emptyInds);
                for kTrial=1:numel(emptyInds)
                    trialRates{emptyInds(kTrial)}=rhat(trialStarts(kTrial):trialEnds(kTrial));
                end
            end
            lambda=cell2mat(trialRates);
        end
        
        %% predict rate within window
        function [rt,bins]=cvPredictAndBinRate(obj, trial, aligningField, win)
            % predict spike rate on witheld data and stitch together
            % rt=cvPredictAndBinRate(obj, trial, aligningField, win)
            nTrials=numel(trial);
            bwin=obj.binfun(win);
            nBins=obj.binfun(diff(win));
            inx=bwin(1):bwin(2)-1;
            if nBins>numel(inx)
                inx=bwin(1):bwin(2); % correct for potential rounding error
            end
            rt=nan(nTrials, nBins);
            
            
            for kFold=1:numel(obj.modelfit)
                rhat=obj.predictSpikeRate(obj.modelfit(kFold), trial, 'trialIndices', obj.modelfit(kFold).testIndices);
                eventTimes=find(obj.getBinnedSpikeTrain(trial, aligningField, obj.modelfit(kFold).testIndices));
                % pad with nans if eventTimes exceeds
                goodix= (eventTimes+bwin(2)<=numel(rhat)) & (eventTimes+bwin(1)>0);
                rt(obj.modelfit(kFold).testIndices(goodix),:)=rhat(bsxfun(@plus, eventTimes(goodix), inx));
            end
            rt=rt(obj.dm.trialIndices,:);
            bins=(bwin(1):bwin(2)-1)*obj.binSize;
        end
        
        
        %% saveobj
        function S=saveobj(obj)
            %             d=obj;
            %             d=funh2struct(obj, false);
            %             S=saveobj@neuroGLM(obj);
            %             S.description=d.description;
            %             S.model=d.model;
            %             S.modelfit=d.modelfit;
            %             S.designOptions=d.designOptions;
            %             S.fitNeuron=d.fitNeuron;
            %             props=properties(obj);
            %             for kProp=1:numel(properties)
            %                 S.(props{kProp})=get(obj, (props{kProp}));
            %             end
            S=obj;
        end
        
        
    end
    
    methods(Static)
        %         function obj=loadobj(s)
        %             d=struct2funh(s,false);
        %             if isstruct(d)
        %                 newObj=glmspike(d.description);
        %                 fields=fieldnames(d);
        %                 for kField=1:numel(fields)
        %                     newObj.(fields{kField})=d.(fields{kField});
        %                 end
        %             else
        %                 newObj=d;
        %             end
        %             obj=newObj;
        %
        %         end
        
        function obj=loadobj(s)
            d=struct2funh(s,false);
            obj=loadobj@neuroGLM(s);
            obj.description=d.description;
            obj.model=d.model;
            obj.modelfit=d.modelfit;
            obj.designOptions=d.designOptions;
            obj.fitNeuron=d.fitNeuron;
        end
        
        %         function obj=loadobj(s)
        %             if isstruct(s)
        %                 newObj=glmspike;
        %                 fields=fieldnames(s);
        %                 for kField=1:numel(fields)
        %                     newObj.(fields{kField})=s.(fields{kField});
        %                 end
        %             else
        %                 newObj=s;
        %             end
        %             obj=struct2funh(newObj,false);
        %         end
    end
    
    
    
    
    
end




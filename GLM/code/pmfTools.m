classdef pmfTools < handle
    % Psychometric Function Fitting Toolbox
    % p = PMFTOOLS(X, Y) - initialize PMFTOOL object
    %   X is a vector of stimulus strengths, Y is a binary vector of
    %   responses. PMFTOOLS automatically bins the stimulus into percentile
    %   bins and calculates the probability of Y given eacy bin
    %
    %
    % p = PMFTOOLS(X,Y,'PARAM1', val1, 'PARAM2', val2, ...)
    %   'nBins'   - number of bins to use when binning data (default = 9)
    %   'binStyle'- how to do stimulus binning
    %               'prctile' (default) nBins using percenttile
    %               'equal'   nBins equally spaced bins
    % p.fit - Maximum Likelihood fit of the data
    % p.fit('PARAM1', val1, ...) allows control of function parameters
    %	'fixedIdx'  - logicals to define which parameters in the fit are
    %                 held fixed and which are free. 
    %                 default is [0 0 1 1] i.e. midpoint & slope are free.
    %                 Fixed parameters get the value set in 'fixedValue':
    %   'fixedValue'- The fixedValue that fixed parameters recieve. 
    %                 default is [0 0] i.e. 0 minLapse and maxLapse
    %
    % p.xValidate - cross validation on withheld data. results stored in
    %   sub-structure 'xvStruct'
    % p.xValidate('PARAM1', val1, ...) allows control of function parameters 
    %	same parameters as p.fit, and:
    %   'kFolds'  - number of folds for cross validation.
    %
    % p.plot - plots the data and the fit
    % p.plot('PARAM1', val1, ...)
    %   'Color'        - Color of the data and fit
    %   'PlotFitBool'  - include the fit in the plot (default = true)
    %   'Marker'       - STRING to set the plot style (default = 'o')
    %
    % --bonus features--
    %   p.visualizeFit    - plots distribution of each parameter
    %   p.test_nTrials    - assesses logL as a function of nTrials
    
    properties
        nTrials     % number of trials
        stim        % input stimulus values
        resp        % input response values
        binId       % id of the bin that each trial belongs to
        x           % x values are binned/discrete stimuli, bootstrapped
        y           % y values are binned/discrete responses, bootstrapped
        err         % [2,1] up & down 1 SEM error, botstrapped
        fitFun      % the fitting function
        nFreeParams % number of free paramteres in the function
        thetaStr    % theta names in string form
        fixedIdx    % the parameters which are fixed
        fixedValue  % the value of the parameters that were fixed.
        theta       % median bootstrapped MLE theta
        thetaErr    % SEM estimated from hessian
        thetaErr95  % the 95% CI
        logL        % median logL
        xFit        % x vlaues for MLE fit
        yFit        % y vlaues for MLE fit
        threshLevel % y values for which you want to compute corresponding x
        threshValue % x values that match threshLevel
        hD          % figure handle for data
        hF          % figure handle for fit
        nDatasets   % optional. if multiple datasets have been combined, this is were to keep track of that.
        testZone    % struct for testing output
        xvStruct  % struct for the xValidate output
    end
    
    
      
            
            
    %%
    methods
        
        %%
        function self = pmfTools(stim, resp, varargin)
            %   stim    - vector of stimulus values
            %   resp    - vector of responses
            %
            % function can handle discrete or continuous stimulus values.
            % if continuous, function bins them into 9 bins.
            
            if nargin < 2
                return;
            end
            
            p = inputParser();
            p.addOptional('nBins', 9);
            p.addOptional('binStyle', 'prctile')
            p.addOptional('binEdges', [])
            p.parse(varargin{:})
            
            stimList    = unique(stim);     % list of stimuli
            nBins       = p.Results.nBins;  % number of bins you wish to bin stim into, if it not discrete
            self.nTrials= size(stim,1);
            self.binId  = nan(self.nTrials,1);
            % If the stimuli are continuous, bin into nBins:
            if numel(stimList) > nBins
                switch p.Results.binStyle
                    case 'equal'
                        [~,binCenters] = hist(stim, nBins);
                        binEdges   = [binCenters(1:2) - diff(binCenters(1:2))/2 binCenters(2:end) + diff(binCenters)/2];
                    case 'prctile'
                        binEdges    = prctile(stim, 0:100/nBins:100);
                        binEdges(1) = binEdges(1) - 0.001; % hack to include 1st element in indices, see iS loop below.
                end
                
                if ~isempty(p.Results.binEdges)
                    binEdges=p.Results.binEdges;
                end
                
                % rounding:
                if binEdges(1) > 5
                    binCenters  = round(binEdges(1:end-1) + diff(binEdges)/2);
                else
                    binCenters  = roundn(binEdges(1:end-1) + diff(binEdges)/2, 2);
                end
                % redefine stimList as the binned stimulus:
                stimList    = binCenters;
                nStim       = numel(stimList);
                pRespMean   = nan(nStim,1);
                pRespErr    = nan(2,nStim);
                % get mean and err of bootstrapped pResp & real:
                for iS = 1:nStim
                    idx = stim > binEdges(iS) & stim <= binEdges(iS+1);
                    self.binId(idx)=iS;
                    pRespMean(iS)=mean(resp(idx));
%                     pRespErr(:,iS)=std(resp(idx))/sqrt(sum(idx));
                    [pRespMean(iS), pRespErr(:,iS)] = bootMeanAndErr(resp(idx));
                end
            else
                % if the stimulus is already discretized:
                nStim       = numel(stimList);
                pRespMean   = nan(nStim,1);
                pRespErr    = nan(2,nStim);
                for iS = 1:nStim
                    idx         = stim==stimList(iS);
                    self.binId(idx)=iS;
                    pRespMean(iS)=mean(resp(idx));
%                     pRespErr(:,iS)=std(resp(idx))/sqrt(sum(idx));
                    [~, pRespErr(:,iS)] = bootMeanAndErr(resp(idx), 2.5, 95);
                end
            end
            % pack up:
            self.stim       = stim;
            self.resp       = resp;
            self.x          = stimList;
            self.y          = pRespMean;
            self.err        = pRespErr;
        end
        
        
        %%
        function self = fit(self, varargin)
            %   self = fit(self, varargin)
            %
            %   fits a logistic function to data via MLE
            %
            % vararg - fixedIdx, fixedValue

            
            p = inputParser();
            % default is to have 2 free parameters, with min/max lapse
            % rate fixed at 0 (i.e. no lapse). For a folded pmf, set the
            % 3rd element (ie 3rd param) in fixedIdx to 1 and send in
            % fixedValue of 0.5.
            p.addOptional('fixedIdx', [0 0 1 1]);
            p.addOptional('fixedValue', [0, 0]);
            p.addOptional('nBoots', 0);
            p.parse(varargin{:})
            
            % number of free parameters:
            self.nFreeParams    = sum(~p.Results.fixedIdx);
            % which parameters will be fixed:
            self.fixedIdx       = p.Results.fixedIdx;
            % their value is determined by fixedValue:
            self.fixedValue     = p.Results.fixedValue;
            % if there are fixed parameters, make sure element numbers match:
            if sum(self.fixedIdx) > 0
                assert(all(sum(self.fixedIdx) == numel(self.fixedValue)), 'ERROR: each fixed parameter must have a fixedValue');
            end
            
            % the psychometric function is defined by the free parameters:
            self.fitFun         = setPsychFun(self.stim, self.fixedIdx, self.fixedValue);
            % names of paramteres:
            self.thetaStr       = {'Inflection', 'Slope', 'Min Lapse', 'Max Lapse'};

            % FIT FUNCTION TO DATA:
            % %%%%%%%%%%%%%%%%%%%%%
            [self.theta, self.logL, self.thetaErr]   = fitMLE(self.stim, self.resp, self.fixedIdx, self.fixedValue);
            self.thetaErr95 = 1.96 .* self.thetaErr; 
            
            % get fit x/y values:
            self.xFit   = linspace(min(self.stim), max(self.stim), 1001);
            fun         = setPsychFun(self.xFit);
            self.yFit   = fun(self.theta);
            
            %%%%%%% BOOTSTRAP OPTION %%%%%%%%
            %%% Not really needed but I'll leave it here as an option, 
            %%% just in case I'll want it in the future...
            if p.Results.nBoots > 1
                nBoots = p.Results.nBoots;
                disp([num2str(nBoots) ' boots will take around ' num2str(nBoots*5e-2) ' seconds to fit...'])
                nVals           = size(self.stim(:),1);
                idxBoot         = RandSample(1:nVals, [nVals, nBoots]);
                thetaDist   = nan(nBoots, 4);
                logLDist    = nan(nBoots, 1);
                for iB = 1:nBoots
                    [thetaDist(iB,:), logLDist(iB)] = fitMLE(self.stim(idxBoot(:,iB)), self.resp(idxBoot(:,iB)), self.fixedIdx, self.fixedValue);
                end
                self.theta      = median(thetaDist);
                self.thetaErr   = mean(abs(bsxfun(@minus, self.theta, prctile(thetaDist, [15.9 84.1])))); % 68 CI
                self.logL       = median(logLDist);
            end
        end
        
        %%
        function self = xValidate(self, varargin)
            %   self = xValidate(self, varargin)
            %
            %   Cross validation of a model on withheld data.
            %   1. Tits a logistic function to 80% of data via MLE
            %   2. Test the likelihood of the model on withheld 20%
            %   3. Repeats for a given number of folds (default 8).
            %
            % vararg - nBoots, nParams, nFolds
            % (100 boots for 2 parms takes approximately 6 seconds)
            
            p = inputParser();
            p.addOptional('kFolds', 8);
            p.parse(varargin{:})
            
            % assert that these data have already been fit. if not, fit:
            if isempty(self.theta)
                self.fit;
            end
                
            % init xvStruct:
            self.xvStruct        = [];
            self.xvStruct.kFolds = p.Results.kFolds;
            % indices for each fold (rows) and train/test sets (columns):
            self.xvStruct.idx    = xvalidationIdx(self.nTrials, self.xvStruct.kFolds, 0, 1);
            
            for k = 1:self.xvStruct.kFolds
                
                % idx for training & test set:
                self.xvStruct.trainIdx(:,k)  = self.xvStruct.idx{k,1}(:);
                self.xvStruct.testIdx(:,k)   = self.xvStruct.idx{k,2}(:);
                % compute pmf on training set and fit:
                self.xvStruct.trainPmf(k)       = pmfTools(self.stim(self.xvStruct.trainIdx(:,k)), self.resp(self.xvStruct.trainIdx(:,k)));
                self.xvStruct.trainPmf(k).fit('fixedIdx', self.fixedIdx, 'fixedValue', self.fixedValue);
                % compute pmf on test set, sans fit:
                self.xvStruct.testPmf(k)       = pmfTools(self.stim(self.xvStruct.testIdx(:,k)), self.resp(self.xvStruct.testIdx(:,k)));
                
                % test training set model likelihood on test data:
                fun = setPsychFun(self.xvStruct.testPmf(k).stim, self.fixedIdx, self.fixedValue);
                y = fun(self.fitFun(self.xvStruct.trainPmf(k).theta));
                self.xvStruct.testL(k) = logLikelihood(self.xvStruct.testPmf(k).resp, y);

            end
        end

        %%
        function self = plot(self, varargin)
            %   self = plot(self, varargin)
            % 
            % plots the binned data along with the fit (if exists)
            % returns handle to data and figure @ self.hF & self.hD
            %
            % varargin - standard plot vars and a few extra: 
            %            PlotFitBool - boolean (0/1) to plot the fit
            %            addLegend - boolean (0/1) to add a legend
            
            
            if isempty(self.x) || isempty(self.y)
                error('ERROR. x & y values are empty. cannot plot data');
            end
            
            figure(gcf); gca;
            fh=plot(nan,nan);
            clr=get(fh, 'Color');
            delete(fh); 
            
            p = inputParser();
            p.addOptional('Color', clr);
            p.addOptional('PlotDataBool', true);
            p.addOptional('PlotFitBool', true);
            p.addOptional('addLegend', false);
            p.addOptional('xLabel', 'Stimulus Value');
            p.addOptional('yLabel', 'Right Choices');
            p.addOptional('fontSize', 10)
            p.addOptional('Marker', '.')
            p.addOptional('MarkerSize', 10)
            p.addOptional('LineWidth', 1)
            p.addOptional('fitLineWidth', 1.5)
            p.addOptional('LineStyle', '-')
            p.parse(varargin{:})
            
            clr          = p.Results.Color;
            marker       = p.Results.Marker;
            markerSize   = p.Results.MarkerSize;
            lineWidth    = p.Results.LineWidth;
            fitLineWidth = p.Results.fitLineWidth;
            lineStyle    = p.Results.LineStyle;
            fontSize     = p.Results.fontSize;
            
            if p.Results.PlotDataBool
                hD = errorbarFancy(self.x, self.y, self.err(1,:), self.err(2,:), ...
                    'Marker', marker, 'MarkerSize', markerSize, 'MarkerEdgeColor',clr, 'MarkerFaceColor', clr, 'Color', clr, 'LineWidth', lineWidth);
            else
                hD  = [];
            end
            hold on;
            
            
            % plot fit, if available:
            if ~isempty(self.xFit) && ~isempty(self.yFit) && p.Results.PlotFitBool
                hF = plot(self.xFit, self.yFit, 'Color', clr, 'LineWidth', fitLineWidth, 'LineStyle', lineStyle);
            else
                hF  = [];
            end
            
            if p.Results.addLegend
                legend([hD hF], {'data', 'fit'}, 'Location', 'NorthWest');
            end
                
            ylim([0 1])
            xlabel(p.Results.xLabel, 'fontSize', fontSize)
            ylabel(p.Results.yLabel, 'fontSize', fontSize)
            
            self.hD = hD; % handle for data
            self.hF = hF; % handle for fit
        end
        
        
        
        %% Extra Functions
        
        %%
        function self = getThreshold(self, threshLevel)
            % get X value of PMF at level (numerical)
            % 
            % input:
            %   threshLevel - y values to obtain x-values (.82 default)
            % output:
            %   threshValue - x values corresponding to y. 
            %   * returns NaN if no value is foudn to match level (e.g. if
            %   the function never crosses the input threshLevel)
            
            if ~exist('threshLevel', 'var')
                threshLevel = 0.82;
            end
            if isempty(self.xFit)
                self.fit;
            end
            
            self.threshLevel = threshLevel;
            for iT = 1:numel(threshLevel)
                % to avoid the interp1 error ("strictly monotically
                % increasing...") I use this hafacilitatory hack:
                ptrStart = find(self.yFit < .95*threshLevel(iT), 1, 'last');
                ptrEnd   = find(self.yFit > 1.05*threshLevel(iT), 1, 'first');
                if ~isempty(ptrStart) && ~isempty(ptrEnd)
                    self.threshValue(iT) = interp1(self.yFit(ptrStart:ptrEnd), self.xFit(ptrStart:ptrEnd), threshLevel(iT));
                else
                    self.threshValue(iT) = nan;
                end
            end

        end
       
        %%
        function self = test_nTrials(self)
            % test how the number of trials affect your fitting
            
            % generate responses from a pmf with known parameters.
            % simulate 5e3 trials from the model.
            % fit MLE to a subset of trials (variable).
            % use fit parms to assess loglikelihood of test set
            
            figure, clf;
            
            % prealloc gen stuct:
            gen         = struct('nTrials', [], 'theta', [], 'stim', [], 'resp', [], 'logL', [], 'xFit', [], 'yFit', []);
            
            % generating model with known parms:
            theta       = [.4, 20 0.05 0.2];
            % simulate responses:
            gen.nTrials     = 5e3;
            gen.stim        = randsample(linspace(0,1,1001)', gen.nTrials, true);
            gen.resp        = rand(gen.nTrials,1) < funpmf('logistic', gen.stim, theta);
            gen             = pmfTools(gen.stim, gen.resp);
            if numel(theta) == 2
                gen = gen.fit;
            elseif numel(theta) == 4
                gen = gen.fit('nParams', 4);
            end
            
            %% define a test-set using the same generating model:
            testSet.nTrials     = 1e3;
            testSet.stim        = randsample(linspace(0,1,1001)', testSet.nTrials, true);
            testSet.resp        = rand(testSet.nTrials,1) < funpmf('logistic', testSet.stim, gen.theta);
            testSet             = pmfTools(testSet.stim, testSet.resp);
            
            %% find best fitting parms to variable nTrials:
            nTrialList  = [100 200 300 500 750];
            nLoops      = numel(nTrialList);
            clear sim;
            for iN = 1:nLoops
                nTrials = nTrialList(iN);
                disp(['fitting ' num2str(nTrials) ' trials'])
                
                % run pmfTools on a subset (nTirals) of the training-data:
                sim(iN) = pmfTools(gen.stim(1:nTrials), gen.resp(1:nTrials));
                % fit the simulated data:
                if numel(gen.theta) == 2
                    sim(iN) = sim(iN).fit;
                elseif numel(gen.theta) == 4
                    sim(iN) = sim(iN).fit('nParams', 4);
                end
                % assess likelihood of the model on the testSet data:
                sim(iN).logL = logLikelihood(testSet.stim, testSet.resp, sim(iN).theta);
                
                %% plot
                % plot simulated data and the MLE fit:
                subplot(2, nLoops, iN); hold on;
                title(sim(iN).nTrials)
                clr = [0 0 0];
                hD = errorbar(testSet.x, testSet.y, testSet.err(1,:), testSet.err(2,:), ...
                    'o', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', clr, 'Color', 'k', 'LineWidth', 1.5);
                % plot the fit to the data:
                hF = plot(sim(iN).xFit, sim(iN).yFit, 'Color', clr, 'LineWidth', 3);
                % plot the generating model:
                hS = plot(gen.xFit, gen.yFit, 'r', 'lineWidth', 1.5);
                
                if iN==1
                    legend([hD, hF, hS], {'test-data', 'fit', 'genModel'}, 'Location', 'NorthWest');
                end
                xlabel('Stimulus Value')
                ylabel('Right Choices')
                ylim([0 1])
                hold off
                hold on;
            end
            % Present logL for fits as a function of nTrials:
            subplot(2, 1, 2); hold on;
            ci      = [sim(:).logLCI];
            ciL     = ci(1:2:end);
            ciU     = ci(2:2:end);
            plot(1:nLoops, [sim(:).logL], '-o', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', clr, 'Color', 'k', 'LineWidth', 1.5);
            title('log likelihood')
            xlim([0.5 nLoops+0.5])
            set(gca,'XTick', 1:nLoops)
            set(gca,'XTickLabel', nTrialList)
            xlabel('nTrials')
            ylabel('log likelihood')
            
            % packup:
            self.testZone.nTrials.gen       = gen;
            self.testZone.nTrials.sim       = sim;
            self.testZone.nTrials.testSet   = testSet;
        end
        
    end
    
    
    methods(Access = private)
        
        
    end
    
end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% Support functions
function y = roundn(x,d)

if nargin == 1
    y = round(x);
else
    y = round(10^d*x)/10^d;
end
end

%% 
function logL = logLikelihood(resp, y)
%   logL = logLikelihood(resp, y)
%
% compute log likelihood of modeled probility 'y' to the responses 'resp'
% INPUT:
%   resp    - 0/1 responses
%   y       - the modeled probability of observing the response
% OUTPUT:
%   logL    - the log likelihood.

% 20151123 - lnk rewrote it

y = y.*0.99999+1e-5; % correction for guessing
logL   = sum(resp.*log(y) + (1-resp).*log(1-y));
end

%%
function psychFun = setPsychFun(stim, fixedIdx, fixedValue)
%   psychFun = setPsychFun(stim, fixedIdx, fixedValue)
%
% the function defines a logistic function of the form:
%   y = gamma + (1-gamma-delta) .* (1./(1+exp(-beta(x-alpha))))
% where each paramter can be set to either free (default) or fixed to
% a value, given by 'fixedValue'. Indices of fixed paramters given by 
% 'fixedIdx'.
%
% example 1:
%   for a folded pmf you might want to set the minLapse to 0.5.
%   Do this by setting fixedIdx=[0 0 1 0] & fixedValue= 0.5  
% exmple 2:
%   for unfolded with fixed lapse rates of 0, set:
%   fixedIdx = [0 0 1 1] & fixedValue = [0 0];
%
% Input:
%   stim        - vector of stimulus values
%   fixedIdx    - vecotr of 0/1 to specify paramters that will remain fixed
%   fixedValue  - vector of parameter values for fixed parameters.
% Output:
%   psychFun    - a function handle to thepsychometric function with your
%                 specified paramter settins.

if exist('fixedIdx', 'var') && ~exist('fixedValue', 'var')
    error('ERROR-- if you fixed some parameters, you must specify their value in fixedValue')
end

if ~exist('fixedIdx', 'var')
    fixedIdx = [0 0 0 0];
end

% all are free:
if all(fixedIdx == [0 0 0 0])
        psychFun    = @(theta) theta(3) + (1-theta(3)-theta(4)) .* (1 ./ (1+exp(-theta(2).*(stim-theta(1)))));
% only one is fixed:        
elseif all(fixedIdx == [0 0 0 1])
        psychFun = @(theta) theta(3) + (1-theta(3)-fixedValue(1)) .* (1 ./ (1+exp(-theta(2).*(stim-theta(1)))));
elseif all(fixedIdx == [0 0 1 0])        
        psychFun = @(theta) fixedValue(1) + (1-fixedValue(1)-theta(3)) .* (1 ./ (1+exp(-theta(2).*(stim-theta(1)))));
elseif all(fixedIdx == [0 1 0 0])
        psychFun = @(theta) theta(2) + (1-theta(2)-theta(3)) .* (1 ./ (1+exp(-fixedValue(1).*(stim-theta(1)))));
elseif all(fixedIdx == [1 0 0 0 ])
        psychFun = @(theta) theta(2) + (1-theta(2)-theta(3)) .* (1 ./ (1+exp(-theta(1).*(stim-fixedValue(1)))));
% two are fixed:
elseif all(fixedIdx == [0 0 1 1])
        psychFun = @(theta) fixedValue(1) + (1-fixedValue(1)-fixedValue(2)) .* (1 ./ (1+exp(-theta(2).*(stim-theta(1)))));
elseif all(fixedIdx == [0 1 0 1])        
        psychFun = @(theta) theta(2) + (1-theta(2)-fixedValue(2)) .* (1 ./ (1+exp(-fixedValue(1).*(stim-theta(1)))));
elseif all(fixedIdx == [1 0 0 1])
        psychFun = @(theta) theta(2) + (1-theta(2)-fixedValue(2)) .* (1 ./ (1+exp(-theta(1).*(stim-fixedValue(1)))));
elseif all(fixedIdx == [0 1 1 0])
        psychFun = @(theta) fixedValue(2) + (1-fixedValue(2)-theta(2)) .* (1 ./ (1+exp(-fixedValue(1).*(stim-theta(1)))));
elseif all(fixedIdx == [1 0 1 0])
        psychFun = @(theta) fixedValue(2) + (1-fixedValue(2)-theta(2)) .* (1 ./ (1+exp(-theta(1).*(stim-fixedValue(1)))));
elseif all(fixedIdx == [1 1 0 0])
        psychFun = @(theta) theta(1) + (1-theta(1)-theta(2)) .* (1 ./ (1+exp(-fixedValue(2).*(stim-fixedValue(1)))));
% three are fixed:
elseif all(fixedIdx == [0 1 1 1])
        psychFun = @(theta) fixedValue(2) + (1-fixedValue(2)-fixedValue(3)) .* (1 ./ (1+exp(-fixedValue(1).*(stim-theta(1)))));
elseif all(fixedIdx == [1 0 1 1])        
        psychFun = @(theta) fixedValue(2) + (1-fixedValue(2)-fixedValue(3)) .* (1 ./ (1+exp(-theta(1).*(stim-fixedValue(1)))));
elseif all(fixedIdx == [1 1 0 1])
        psychFun = @(theta) theta(1) + (1-theta(1)-fixedValue(3)) .* (1 ./ (1+exp(-fixedValue(2).*(stim-fixedValue(1)))));
elseif all(fixedIdx == [1 1 1 0])
        psychFun = @(theta) fixedValue(3) + (1-fixedValue(3)-theta(1)) .* (1 ./ (1+exp(-fixedValue(2).*(stim-fixedValue(1)))));
else
        error('unsupported')
end

end


    

%%
function [theta, logL, SEM] = fitMLE(stim, resp, fixedIdx, fixedValue, thetaConLo, thetaConHi)
%   [theta, logL] = fitMLE(stim, resp, freeIdx, fixedValue, thetaConLo, thetaConHi)
%
% Fits a Maximum Likelihood Estimate (MLE) of the to input data using a
% logistic function.
%
% INPUTS:
%   cond        - vecotr of conditions (e.g. stim intensity)
%   resp        - vecotr of responses (0s and 1s)
%   fixedIdx    - 0/1 determines which parameters are fixed and which are 
%                 free. Fixed params take their value from theta0.            
%   fixedValue  - an array of size sum(fixedIdx) that has a value for each
%                 parameter that is fixed.
%   thetaConLo  - low constraint for the theta
%   thetaConHi  - high constraint for theta
% OUTPUTS:
%   theta       - theta that minimizes the negative log likelihood
%   logL        - log likelihod produced by theta
%   SEM         - The std of the posterior estimate of the parameter, equal
%                 to 68.2% CI (i.e. SEM) derived from the numerically 
%                 computed hessian (muliply by 1.96 if you want the 95%)

% 2014/06 lnk wrote it.
% 201504 lnk optimized it.
% 201511 lnk added hessian output
% 201511 lnk added fixedIdx to control which parameters are fixed.

% More info:
% this function relies on setPsychFun. setPsychFun defines the logisitc 
% function that will be fit. Its paramteres can be controlled by fixedIdx.
% Thus, the objective function already has the free paramteres to minize
% over.
% for any fixed paramter (i.e. a 'true' value in fixedIdx) the value the
% parameter will take is the same as the correspoiding theta0 parameter.


if ~exist('fixedIdx', 'var')
    fixedIdx = [0 0 0 0];
end
fixedIdx    = logical(fixedIdx);
nParams     = numel(fixedIdx);

if ~exist('thetaConLo', 'var') || ~exist('thetaConHi', 'var')
    % default theta constraints (talilored for logistic fits. wbl would be differnt):
    thetaConLo = [min(stim),  1 / range(stim),   0,    0];
    thetaConHi = [max(stim),  100 / range(stim), .6,     .6];
end

% First, define appropriate function given that number of free parameters.
% Next,  set objective function to be minimized as the NLL:
psychFun    = setPsychFun(stim, fixedIdx, fixedValue);
objFun      = @(theta) (-1)*logLikelihood(resp, psychFun(theta));

opts        = optimset('display', 'none', 'algorithm', 'interior-point');

%% guess starting parameters:
% build a grid of all non-fixed paramter permutations, and assess logL.
% Choose combination with highest logL as 'thetaGuess' for minimization.
pGrid = cell(nParams,1);
id = 1;
for ii = 1:nParams
    if fixedIdx(ii)
        pGrid{ii} = fixedValue(id);
        id = id+1;
    else
        pGrid{ii} = linspace(thetaConLo(ii), thetaConHi(ii), 10);
    end
end
        
% create grid for all possible combinations:
[p1, p2, p3, p4]    = ndgrid(pGrid{1}, pGrid{2}, pGrid{3}, pGrid{4});
phiGrid             = [p1(:), p2(:), p3(:), p4(:)];
% removing instances were p3+p4 is more than 1:
phiGrid = phiGrid(p3(:)+p4(:) <= 1,:);

% loop through each combo and get the logLikelihood:
nPhi        = size(phiGrid,1);
logL        = nan(nPhi, 1);
for ii = 1:nPhi
    logL(ii) = logLikelihood(resp, psychFun(phiGrid(ii,:)));
end
ptr         = find(logL==max(logL),1);
thetaGuess  = phiGrid(ptr,:);
logLGuess   = logLikelihood(resp, psychFun(thetaGuess));

% minimize the objective function starting from the best guess 'thetaGuess'
% only for un-fixed paramteres:
[thetaMLE, ~, exitflag] = fmincon(objFun, thetaGuess(~fixedIdx), [],[],[],[], thetaConLo(~fixedIdx), thetaConHi(~fixedIdx), [], opts);
% compute the hessian numerically and estimate the error:
H       = compHess(objFun, thetaMLE', 0.001);
tmpSEM  = sqrt(diag(inv(H)));

% and then re-introduce fixed paramteres into correct place:
theta   = nan(1, numel(fixedIdx));
SEM     = nan(1, numel(fixedIdx));
iMLE    = 1;
iFixed  = 1;
if sum(fixedIdx) > 0
    for ii = 1:nParams
        if fixedIdx(ii)
            theta(ii)   = fixedValue(iFixed);
            SEM(ii)     = fixedValue(iFixed);
            iFixed = iFixed+1;
        else
            theta(ii)   = thetaMLE(iMLE);
            SEM(ii)     = tmpSEM(iMLE);
            iMLE = iMLE+1;
        end
    end
else
    theta   = thetaMLE;
    SEM     = tmpSEM;
end



% and evaluate logL:
logL        = logLikelihood(resp, psychFun(theta));

if logLGuess > logL
    warning('logMLE failed to maximize likelihood beyond that obtained by the grid. USER-- solve this')
    keyboard;
end
if exitflag <= 0
    theta = nan(size(theta));
    warning('Couldn''t find a minima, exiting without returning best fit parameters');
    keyboard;
end

end

%%
function [bootMean, bootErr, bootDstMeans] = bootMeanAndErr(dst, lPrc, uPrc, nBoots)
%   [bootMean, bootErr, bootDstMeans] = bootMeanAndErr(dst, lprc, uprc)
%
% bootstrasp a given dstribution to compute a mean and error bar.
% Errorbar is set by lPrc & uPrc, default is 1 SEM, ie 68.2% confidence
% interval.
%
% Inputs:
%   dst         - the dstribution of numbers.
%   [lprc]      - lower percentile (default is 15.9)
%   [uprc]      - upper percentile (default is 84.1)
%   [nBoots]    - number of boots (default is 1000)
% Outputs:
%   bootMean        - mean of means of bootstrapped dstributions
%   bootErr         - [2,1] lower and upper values for error bars
%   bootDstMeans   - [nBoots,1] the means for all bootstrapped dsts.
%
% See also bootMeanAndCI

% 201406 lnk wrote it

%%

if nargin < 4
    nBoots = 1e3;
    if nargin < 3
        lPrc = 15.9;
        uPrc = 84.1;
    end
end

%%

if ~isempty(dst)
    nValues = numel(dst);
    
    bootDst        = reshape(randsample(dst, nValues*nBoots,true),[nValues, nBoots]);
    bootDstMeans   = nanmean(bootDst,1);
    
    bootMean        = nanmean(bootDstMeans);
    bootCI          = prctile(bootDstMeans, [lPrc, uPrc]);
    bootErr         = abs(bootMean - bootCI);
else
    bootMean    = nan;
    bootErr     = [nan,nan];
end

end


%%
function h = errorbarFancy(x,y,l,varargin)
%   h = errorbarFancy(x,y,l,u,varargin)
%
%   Error bar plot without the annoying horizontal lines.
%   Based on drawing a line, and then plotting.
%   takes in same arguments as the plot function, see LineSpec.

if nargin > 3 && isnumeric(varargin{1})
    u=varargin{1};
    if numel(varargin)>1
        varargin = varargin(2:end);
    else
        varargin={};
    end
else
    u=l;
end

% if ~exist('u', 'var') && exist('l', 'var')
%     u = l; % symmetrical errorbars
% end

% hack to find a lineWidth, if inputted:
if any(strcmp(varargin,'LineWidth'))
    idx = find(strcmp(varargin,'LineWidth'));
    lw  = varargin{idx+1};
else
    lw = 1; % default;
end

% hack to find a Color, if inputted:
if any(strcmp(varargin,'Color'))
    idx = find(strcmp(varargin,'Color'));
    clr  = varargin{idx+1};
else
    fh=plot(nan,nan);
    clr=get(fh, 'Color');
    delete(fh); 
    varargin={varargin{:}, 'Color', clr};
    %     clr = [0 0 0]; % default;
end

gca; hold on;
line([x(:) x(:)]', [y(:)-l(:) y(:)+u(:)]', 'Color', clr, 'LineWidth', lw);
h = plot(x,y, 'Marker','o', 'LineStyle', 'None', varargin{:});


end

%%
function xidxs = xvalidationIdx(nTotalSamples, kxValidation, isRandomized, isPartition)
% xidxs = xvalidationIdx(nTotalSamples, kxValidation, isRandomized, isPartition)
% Get k-fold cross-validation indices.
%
% Input
%   nTotalSamples: integer number of samples or indices
%   kxValidation: k for the k-cross validation (if 1, training = test)
%   isRandomized: (opt) sequencial or randomized (default) indices?
%   isPartition: (opt) drop extra points at the end to make equal size partitions
% Output
%   xidxs: {kxValidation x 2} training, test indices
%
% Caution: test set is not guarantteed 
%          to be a partition, UNLESS mod(nTotalSamples, kxValidation) == 0
%          In such case, the test set is a partition of nTotalSamples.
%	   Use isPartition to enforce partitioning.
%
% $Id$
% Copyright 2011 Memming. All rights reserved.

if nargin < 3
    isRandomized = true;
end

if nargin < 4
    isPartition = false;
end

if isPartition
    nTotalSamples = nTotalSamples - rem(nTotalSamples, kxValidation);
end

if rem(nTotalSamples,1) ~= 0 || rem(kxValidation,1) ~= 0
    error('xvalidationIdx', 'arguments should be both integers');
end

if ~numel(nTotalSamples) == 1
    ridx = nTotalSamples;
    nTotalSamples = numel(ridx);
else
    ridx = 1:nTotalSamples;
end

% size of the test set
m = ceil(nTotalSamples / kxValidation);
if m < 1
    error('xvalidationIdx:insufficient_samples', 'Not enough samples (%d) to create %d-fold cross-validation', nTotalSamples, kxValidation);
end

if isRandomized
    sidx = randperm(nTotalSamples);
    ridx = ridx(sidx);
end

if kxValidation == 1
    % No cross-validation. All samples are training, and all samples are test.
    xidxs{1,1} = ridx;
    xidxs{1,2} = ridx;
    return
end

startIdx = ceil(linspace(1, nTotalSamples, kxValidation+1));
for k = 1:kxValidation
    testSet = startIdx(k) + (1:m) - 1;
    trainSet = setdiff(1:nTotalSamples, testSet);
    xidxs{k,1} = ridx(trainSet);
    xidxs{k,2} = ridx(testSet);
end

end

%%  
function [H, g] = compHess(fun, x0, dx, varargin)
% [H, g] = compHess(fun, x0, dx, varargin)
% Numerically computes the Hessian of a function fun around point x0
% expects fun to have sytax:  y = fun(x, varargin);
%
% Input:
%   fun: @(x) function handle of a real valued function that takes column vector
%   x0: (n x 1) point at which Hessian and gradient are estimated
%   dx: (1) or (n x 1) step size for finite difference
%   extra arguments are passed to the fun
%
% Output:
%   H: Hessian estimate
%   g: gradient estiamte

n = numel(x0);
H = zeros(n,n);
g = zeros(n,1);
f0 = feval(fun, x0, varargin{:});

% input check
if ~all(isfloat(dx) & isfinite(dx))
    error('dx must be finite and float');
end
if any(dx <= 0)
    error('dx must be strictly positive');
end

if isscalar(dx)
    vdx = dx*ones(n,1);
elseif numel(dx) == n
    vdx = dx(:);
else
    error('vector dx must be the same size as x0');
end
A = diag(vdx/2);

for j = 1:n  % compute diagonal terms
    %central differences
    f1 = feval(fun, x0+2*A(:,j),varargin{:});
    f2 = feval(fun, x0-2*A(:,j),varargin{:});
    H(j,j) = f1+f2-2*f0;
    g(j) = (f1-f2)/2;
    % forward differences
%     f1 = feval(fun, x0+2*A(:,j),varargin{:});
%     f2 = feval(fun, x0+A(:,j),varargin{:});
%     fx = feval(fun, x0,varargin{:});
%     H(j,j) = f1-2*f2+fx;
%     g(j)   = f2-fx;
end

for j = 1:n-1       % compute cross terms
    for i = j+1:n
        %central differences
        f11 = feval(fun, x0+A(:,j)+A(:,i),varargin{:});
        f22 = feval(fun, x0-A(:,j)-A(:,i),varargin{:});
        f12 = feval(fun, x0+A(:,j)-A(:,i),varargin{:});
        f21 = feval(fun, x0-A(:,j)+A(:,i),varargin{:});
        H(j,i) = f11+f22-f12-f21;
        % forward differences
%         fx = feval(fun, x0,varargin{:});
%         f1 = feval(fun, x0+A(:,j)+A(:,i),varargin{:});
%         f12 = feval(fun, x0+A(:,j),varargin{:});
%         f21 = feval(fun, x0+A(:,i),varargin{:});
%         H(j,i) = f1-f12-f21-fx;

        H(i,j) = H(j,i);
    end
end

H = H./(vdx * vdx');
g = g./vdx;

end

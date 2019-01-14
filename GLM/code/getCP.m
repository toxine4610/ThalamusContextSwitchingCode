function CP = getCP(stim, spikeTimes, trialIndices, targPref, verbose)
% CP = getCP(stim, spikeTimes, trialIndices, targPref, verbose)

if ~exist('verbose', 'var')
    verbose = false;
end

if ~exist('targPref', 'var')
   targPref = 1; 
end

if ~exist('trialIndices', 'var')
    trialIndices = stim.goodtrial & stim.dirprob == 0;
end

opts = struct('aligningWindow', [-1 2], ...
    'binSize', .01, ...
    'smoothingKern', ones(10, 1)/10);

CP = struct();
CP.trialIndices = trialIndices;

% centeringField = stim.timing.motionOnset(:,1) + stim.timing.plxstart;
centeringField = stim.timing.motionon(:,1) + stim.timing.plxstart;

[~,~,binCenters, ~, trialSpikeCount] = eventPsth(spikeTimes, centeringField(trialIndices), opts.aligningWindow, opts.binSize, opts.smoothingKern);

X = sum(stim.pulses,3); 
targPlus = mode(stim.targchosen(sign(sum(X,2)) == stim.targchosen));

CP.idx = trialIndices;

if all(stim.frozentrials(trialIndices))
    S = X(trialIndices,:);
    R = stim.targchosen(trialIndices) == targPref;
    [uniqueStimSequences, ~, IAB] = unique(S,'rows');
    nConds = size(uniqueStimSequences,1);
    nBins = numel(binCenters);
    
    m = zeros(nBins, nConds); s = m;
    ms = zeros(stim.nPulses, nConds); sm = ms;
    n = zeros(nConds,1);
    for kFrozen = 1:nConds
        n(kFrozen) = sum(IAB==kFrozen);
        [m(:,kFrozen), s(:, kFrozen)] = rcg.choiceProbabilityCalculate(trialSpikeCount(IAB==kFrozen,:), R(IAB==kFrozen));
        [ms(:,kFrozen), sm(:, kFrozen)] = rcg.choiceProbabilityCalculate(S(IAB==kFrozen,:), R(IAB==kFrozen));
    end
    
    CP.cpBinCenters 		= binCenters;
    CP.cpMean  			= (m*n)/sum(n);
    CP.cpError 			= (s*n)/sum(n);
    CP.pulseBinCenters 	= (0:stim.nPulses-1) * stim.nFrames/stim.nPulses/stim.frate;
    CP.pulseAUCMean 		= (ms*n)/sum(n);
    CP.pulseAUCErr  		= (sm*n)/sum(n);
    CP.nTrials 			= sum(n);
else
    
    [ms, sm] = rcg.choiceProbabilityCalculate(X(CP.idx,:), stim.targchosen(CP.idx)==targPlus);
    [m, s]   = rcg.choiceProbabilityCalculate(trialSpikeCount, stim.targchosen(CP.idx)==targPref);
    
    CP.cpBinCenters 		= binCenters;
    CP.cpMean  			= m;
    CP.cpError 		    = s;
    CP.pulseBinCenters    = (0:stim.nPulses-1) * stim.nFrames/stim.nPulses/stim.frate;
    CP.pulseAUCMean 		= ms;
    CP.pulseAUCErr  		= sm;
    CP.nTrials 			= sum(CP.idx);
    
end

if verbose
    figure(kCondition); clf
    errorbar(CP.pulseBinCenters, CP.pulseAUCMean, CP.pulseAUCErr, 'k'); hold on
    errorbar(CP.cpBinCenters, CP.cpMean, CP.cpError, 'b')
end

end


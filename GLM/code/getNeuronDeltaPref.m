function deltaPref = getNeuronDeltaPref(neuron, stim)
% get delta pref for a neuron on a specific session
% deltaPref = getNeuronDeltaPref(neuron, stim)

assert(~iscell(neuron), 'neuron must be a struct')

prefDir = getNeuronPrefDir(neuron);

deltaPref = nan;

if isnan(prefDir)
    return
end

switch neuron.brainArea
    case 'MT'
        deltaPref = min(getDeltaTuning(stim.theta, prefDir), getDeltaTuning(stim.theta-180, prefDir));
    case 'LIP'
        targ1XY = mode(stim.targ1XY);
        targ2XY = mode(stim.targ2XY);
        targ1Th = cart2pol(targ1XY(1), targ1XY(2))/pi*180;
        targ2Th = cart2pol(targ2XY(1), targ2XY(2))/pi*180;
        
        deltaPref = min(getDeltaTuning(targ1Th, prefDir), getDeltaTuning(targ2Th, prefDir));
end
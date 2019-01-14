function prefDir = getNeuronPrefDir(neuron)
% get a neuron's preferred direction
% prefDir = getNeuronPrefDir(neuron)
% Inputs
%   neuron (struct) - output of make_SingleNeuronStruct
% Outputs
%   prefDir - in degrees

prefDir = nan;

switch neuron.brainArea
    case 'MT'
        hyperflowPrefDir = nan;
        mtrfmapPrefDir   = nan;
        
        if ~isempty(neuron.mtrfmap)
            mtrfmapPrefDir = neuron.mtrfmap.prefDir;
            prefDir=mtrfmapPrefDir;
            return
        end
        
        if ~isempty(neuron.hyperflow) && isfield(neuron.hyperflow, 'prefDir')
            try
                hyperflowPrefDir = neuron.hyperflow.prefDirInRF;
            catch
                hyperflowPrefDir = neuron.hyperflow.prefDir;
            end
            prefDir=hyperflowPrefDir;
            return
        end
        
        if ~isnan(hyperflowPrefDir) && ~isnan(mtrfmapPrefDir)
            assert(getDeltaTuning(mtrfmapPrefDir, hyperflowPrefDir) < 20, 'the two mapping styles give different answers! What is up?')
        end
        
%         prefDir = nanmean([hyperflowPrefDir mtrfmapPrefDir]);
        
        
    case 'LIP'
        
        if ~isempty(neuron.delayedsaccades)
            prefDir = cart2pol(neuron.delayedsaccades.rfXY(1), neuron.delayedsaccades.rfXY(2))*180/pi;
        end
        
    otherwise
        warning('getNeuronPrefDir: no support for areas besides MT and LIP')
end
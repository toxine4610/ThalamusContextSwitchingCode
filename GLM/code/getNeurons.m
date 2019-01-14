function neurons = getNeurons(exname, dataPath)
% GETNEURONS gets spike times and meta data for experiment <exname>
% INPUT: 
% 	exname   <char>  experiment name
%   dataDir  <path>  full path to data
% OUPUT:
% 	neurons <neuro>  array of neuron objects
%
% EXAMPLE CALL
% 	neurons = getNeurons('p20140304', './data');

fl = dir(fullfile(dataPath, 'neurons', [exname '*']));

if isempty(fl)
    neurons = [];
    warning('Neuron data file corresponding to [%s] not found. Check [%s]', exname, fullfile(dataPath, 'neurons'));
    return
end

for k = 1:numel(fl)
    neurons(k) = neuro(fullfile(dataPath, 'neurons', fl(k).name));
end
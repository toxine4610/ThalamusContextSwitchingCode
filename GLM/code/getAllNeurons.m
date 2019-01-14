function nw = getAllNeurons(dataPath, tags)
% nw = getAllNeurons(tags)

if ~exist('tags', 'var')
    tags = {'good', 'simultaneous'};
end
    
experiments = getExperimentsAnd(tags);
nExperiments = numel(experiments);

for kEx = 1:nExperiments
    exname = experiments{kEx};
    if kEx == 1
        nw = getNeurons(exname, dataPath);
        continue
    end
    
    nw = [nw getNeurons(exname, dataPath)];
end

if isempty(nw)
    error('No neurons loaded! Make sure data files are in [%s]', dataPath);
end
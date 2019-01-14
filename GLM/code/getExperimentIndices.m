function [Indices, experiments] = getExperimentIndices(tags)
% get a set of indices for experiments that match <tags>
% [Indices, experiments] = getExperimentIndices(tags)
% example:
%  Indices, experiments = getExperimentIndices{'MT', 'LIP', 'GoodLFP'}
% 	returns a logical index of size [nExperiments x nTags]

% jly wrote it

evalc('[~, experiments] = GaborDotsNotes();');

nExperiments = numel(experiments);

nTags = numel(tags);
Indices = false(nExperiments, nTags);

for kExperiment = 1:nExperiments
    
    metaData = GaborDotsNotes(experiments{kExperiment});
    
    for kTag = 1:nTags
        if isfield(metaData, 'tags')
            Indices(kExperiment, kTag) = any(strcmp(metaData.tags, tags{kTag}));
        end
        
    end
end


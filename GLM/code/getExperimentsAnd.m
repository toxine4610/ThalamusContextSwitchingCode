function experiments = getExperimentsAnd(tags)
% get experiment names based on tags (performs AND operation on all matches)
if ischar(tags)
    tags = {tags};
end

[~,experiments] = GaborDotsNotes();
good = getExperimentIndices(tags);
if size(good,2) > 1
    good = all(good,2);
end
experiments = experiments(good);
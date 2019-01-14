function S = loadUpFits(dataPath, tag, fitDir)
% S = loadUpFits(dataPath, tag)

if nargin < 3
    fitDir = 'main_fits';
end

assert(isdir(dataPath), 'First argument must be the path to data')

nw = getAllNeurons(dataPath);

ixDprime = abs([nw(:).dPrimeSigned]) > .2;
ixLIP    = [nw(:).isLIP];

switch tag
    case 'MT'
        ix = ixDprime & ~ixLIP;
    case 'LIP'
        ix = ixDprime & ixLIP;
    otherwise
        ix = ixDprime;
end

nw(~ix) = [];

clear S
for k = 1:numel(nw)
    fname = fullfile(dataPath, fitDir, [nw(k).getName 'modelComparison.mat']);
    tmp = load(fname);
    S(k) = tmp.model;
end
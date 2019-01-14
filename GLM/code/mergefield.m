function [dv] = mergefield(dv,defaults,nowarning)
% merge one struct into another (overwride when duplicate values)
% [dv] = mergefield(dv,defaults,nowarning)
% merge dv->defaults replacing fields of defaults with values specified in
% dv

fn = fieldnames(defaults);
for ii = 1:length(fn)
    if ~isfield(dv,fn{ii}) || isempty(dv.(fn{ii}))
        dv.(fn{ii}) = defaults.(fn{ii});
    end
end

if nargin < 3 || nowarning == 0
    fn = fieldnames(dv);
    for ii = 1:length(fn)
        if ~isfield(defaults,fn{ii});
            warning('setdefaults:unknownarg','Argument %s is unsupported',fn{ii});
        end
    end
end

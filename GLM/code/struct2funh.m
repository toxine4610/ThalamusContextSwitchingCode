function obj=struct2funh(obj, verbose)
% replace function meta data with appropriate function handle
% obj=struct2funh(obj, verbose, recursive)
% verbose true prints out every function handle that is repaired

if ~exist('verbose', 'var') || isempty(verbose)
    verbose=true;
end

recursive=true;
if isstruct(obj)
    fields=fieldnames(obj);
else
    fields=properties(obj);
end

for kField=1:numel(fields)
    if isnumeric(obj.(fields{kField})) || ischar(obj.(fields{kField}))
        continue
    end
    n=numel(obj.(fields{kField}));
    if n>1
        for k=1:n
            if isa(obj.(fields{kField})(k), 'struct') && isfield(obj.(fields{kField})(k), 'workspace')
                obj.(fields{kField})(k)=reconstructFcn(obj.(fields{kField})(k));
            elseif (isa(obj.(fields{kField})(k), 'handle') || (isa(obj.(fields{kField})(k), 'struct') && ~isfield(obj.(fields{kField})(k), 'workspace'))) && recursive
                obj.(fields{kField})(k)=struct2funh(obj.(fields{kField})(k), verbose);
            end
        end
    else
        if isa(obj.(fields{kField}), 'struct') && isfield(obj.(fields{kField}), 'workspace')
            obj.(fields{kField})=reconstructFcn(obj.(fields{kField}));
            if verbose
                fprintf('Function [%s] reconstructed\n', (fields{kField}))
            end
        elseif (isa(obj.(fields{kField}), 'handle') || (isa(obj.(fields{kField}), 'struct') && ~isfield(obj.(fields{kField}), 'workspace'))) && recursive
            obj.(fields{kField})=struct2funh(obj.(fields{kField}), verbose);
        end
    end
end
if verbose
    fprintf('Done\n')
end

function out = reconstructFcn(s)

for iWks = 1:numel(s.workspace)
    wkspace = s.workspace{iWks};
    varnames = fieldnames(wkspace);
    if ~isempty(varnames)
        for i = 1:numel(varnames)
            try %#ok<TRYNC>
                tmp = wkspace.(varnames{i}); %#ok<NASGU>
                eval([varnames{i} ' = tmp;']);
            end
        end
    end
end

fcn_str = s.function;
fcn_str=fcn_str(strfind(fcn_str, '@'):end);
fcn = eval(fcn_str);
out = fcn;



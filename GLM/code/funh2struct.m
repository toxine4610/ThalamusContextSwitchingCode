function obj=funh2struct(obj, verbose)
% convert all function handles in the hierarchy to structs with meta data
% obj=funh2struct(obj, verbose, recursive)
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
            if (isa(obj.(fields{kField})(k), 'handle') || (isa(obj.(fields{kField})(k), 'struct') && ~isfield(obj.(fields{kField})(k), 'workspace'))) && recursive
                obj.(fields{kField})(k)=funh2struct(obj.(fields{kField})(k), verbose);
            end
        end
    else
       if (isa(obj.(fields{kField}), 'handle') || (isa(obj.(fields{kField}), 'struct') && ~isfield(obj.(fields{kField}), 'workspace'))) && recursive
            obj.(fields{kField})=funh2struct(obj.(fields{kField}), verbose);
        elseif isa(obj.(fields{kField}), 'function_handle')
            f=functions(obj.(fields{kField}));
            if strcmp(f.type, 'anonymous')
                obj.(fields{kField})=f;
                if verbose
                    fprintf('Function handle [%s] meta data included\n', fields{kField})
                end
            end
        end
    end
end
if verbose
    fprintf('Done\n')
end



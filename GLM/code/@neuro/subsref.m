function varargout = subsref(N,S)
            % handle the case of .calls
            if  strcmp(S(1).type, '.')
                if any(strcmp(S(1).subs,methods(N))) || any(strcmp(S(1).subs,properties(N)))
                    % Enable dot notation for some functions
                    [varargout{1:nargout}] = builtin('subsref',N,S);
                else
                    [varargout{1:nargout}] = builtin('subsref',N.nreader,S);
                end
                return
            end
            
            if numel(S)>1 && strcmp(S(1).type, '()')
                n = numel(N(S(1).subs{:}));
                if ischar(S(1).subs{:}) && n > 1% check for colon
                    ctr = 1;
                    for k = 1:n
                        [varargout{ctr}] = N(k).subsref(S(2:end)); %#ok<AGROW>
                        ctr = ctr+1;
                    end
                    if nargout == 1 % overload for arrays of neuros
                        varargout = {[varargout{:}]};
                    end
                elseif n > 1
                    ctr = 1;
                    for k = S(1).subs{:}(:)'
                        [varargout{ctr}] = N(k).subsref(S(2:end)); %#ok<AGROW>
                        ctr = ctr+1;
                    end
                    if nargout == 1 % overload for arrays of neuros
                        varargout = {[varargout{:}]};
                    end
                else
                    for k = S(1).subs{:}(:)'
                        [varargout{1:nargout}] = N(k).subsref(S(2:end));
                    end
                end
            end
            
            % if nothing is passed in enable standard notation
            if numel(S)==1 && strcmp(S(1).type, '()')
                [varargout{1:nargout}] = builtin('subsref',N,S);
            end
            
        end
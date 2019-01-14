% RCGREADER loads stimulus properties from disk when needed
classdef rcgreader
    properties
        exname
        filename
        rdr
    end
    
    methods
        function s = rcgreader(exname, dataPath)
            fname=fullfile(dataPath, 'stim', [exname '_stim.mat']);
            s.exname    = exname;
            s.filename  = fname;
            s.rdr       = reader(fname);
        end
        
        function t = targPlus(s)
            % round sum of pulses to 1 or 2 (target id)
            xsign = sign(sum(sum(s.rdr.pulses,3),2));
            xzero = xsign == 0;

            pos = round(1.5-.5*xsign);
    
            propMatchTarg1Plus = mean(pos==s.rdr.targcorrect | xzero);
            propMatchTarg2Plus = mean(pos~=s.rdr.targcorrect | xzero);

            if propMatchTarg1Plus > propMatchTarg2Plus
                t = 1;
            else
                t = 2;
            end

        end
        
        
        function varargout = subsref(s,S)
            
            switch S(1).type
                case '.'
                   if any(strcmp(S(1).subs,methods(s))) || any(strcmp(S(1).subs,properties(s)))
                        % Enable dot notation for some functions
                        if(nargout==0)
                            builtin('subsref',s,S);
                        else   
                            [varargout{1:nargout}] = builtin('subsref',s,S);%p.(S.subs);  
                        end
                   elseif any(strcmp(S(1).subs,'timing'))
                       tim = s.rdr.timing;
                       
                       if numel(S) == 1
                           varargout{1} = tim;
                           return
                       end
                       
                        if ~ischar(S(2).subs)
                            b = tim(S(2).subs{:});
                            if numel(S) > 2
                                b = [b.(S(3).subs)];
                            end
                            varargout{1} = b;
                            return
                        end
                       
                       nTrials = numel(tim);
                       m = nan(nTrials, 1);
                       for kTrial = 1:nTrials
                           m(kTrial) = numel(tim(kTrial).(S(2).subs));
                       end
                           
                       b = nan(nTrials, max(m));
                       for kTrial = 1:nTrials
                           tmp = tim(kTrial).(S(2).subs);
                           b(kTrial,1:numel(tmp)) = tmp;
                       end
                       
                       if numel(S)>2
                           varargout{1} = b(S(3).subs{:});
                       else
                           varargout{1} = b;
                       end
                       
                   else
                       if(nargout==0)
                           builtin('subsref',s.rdr,S);
                       else   
                           [varargout{1:nargout}] = builtin('subsref',s.rdr,S);%p.(S.subs);  
                       end
                       
                   end
            end
       end
    end
    
end
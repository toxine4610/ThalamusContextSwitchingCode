classdef reader
    properties
        filename
        st
        isFile
        isStruct
    end
    
     methods 

       function r=reader(in)
           if nargin > 0
               if isstruct(in)
                   r.st=in;
                   r.filename=[];
                   r.isFile=false;
                   r.isStruct=true;
               elseif exist(in, 'file')
                   r.filename=in;
                   r.st=[];
                   r.isFile=true;
                   r.isStruct=false;
               end
           else
               r.filename = [];
               r.st = [];
               r.isFile = [];
               r.isStruct=[];
           end
       end
       
       function varargout = subsref(r,S)
            
            switch S(1).type
                case '.'
                   if any(strcmp(S(1).subs,methods(r)))
                        % Enable dot notation for some functions
                        [varargout{1:nargout}] = builtin('subsref',r,S);
                   elseif r.isFile
                       tmp=load(r.filename,S(1).subs);
                       [varargout{1:nargout}] = builtin('subsref',tmp,S);
                   elseif r.isStruct
                       [varargout{1:nargout}] = builtin('subsref',r.st,S);
                   else
                       [varargout{1:nargout}] = nan;
                   end
            end
       end
                       
       
     end
     
end
          
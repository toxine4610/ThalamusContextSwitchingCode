
function DelT = computeTimeStamps(LEDstate, timestamps);

s = strfind(LEDstate, 'autostart');
index = find(~cellfun(@isempty,s));
isrepeat = []; norepeat = []; success = [];

isbrown          = zeros(1,length(index));
isblue           = zeros(1,length(index));
% isrepeat               = zeros(1,length(index));
% % isrepeat_same     = zeros(1,length(index));
% norepeat          = zeros(1,length(index));
%
for i = 1:length(index)
    curr_index = index(i);
    
    if strcmpi( LEDstate{curr_index-2}, 'LEDTG On' ) || strcmpi( LEDstate{curr_index-2}, 'LEDTB On' )
        issound(i) = 1;
        DelT(i)    = str2num( timestamps{ curr_index+1 } ) - str2num( timestamps{ curr_index-2+1 } );
    elseif strcmpi( LEDstate(curr_index-1), 'red should flash...') && strcmpi( LEDstate(curr_index-5), 'LEDTG On')
        issound(i) = 1;
        DelT(i)    = str2num( timestamps{ curr_index+1 } ) - str2num( timestamps{ curr_index-5+1 } );
    elseif strcmpi( LEDstate(curr_index-1), 'red should flash...') && strcmpi( LEDstate(curr_index-5), 'LEDTB On')
        issound(i) = 1;
        DelT(i)    = str2num( timestamps{ curr_index+1 } ) - str2num( timestamps{ curr_index-5+1 } );
    
    elseif strcmpi( LEDstate(curr_index-1), 'red should flash...') && strcmpi( LEDstate(curr_index-4), 'LEDTG On')
        issound(i) = 1;
        DelT(i)    = str2num( timestamps{ curr_index+1 } ) - str2num( timestamps{ curr_index-4+1 } );
    elseif strcmpi( LEDstate(curr_index-1), 'red should flash...') && strcmpi( LEDstate(curr_index-4), 'LEDTB On')
        issound(i) = 1;
        DelT(i)    = str2num( timestamps{ curr_index+1 } ) - str2num( timestamps{ curr_index-4+1 } );
        
    elseif strcmpi( LEDstate(curr_index-1), 'red should flash...') && strcmpi( LEDstate(curr_index-3), 'LEDTG On')
        issound(i) = 1;
        DelT(i)    = str2num( timestamps{ curr_index+1 } ) - str2num( timestamps{ curr_index-3+1 } );
    elseif strcmpi( LEDstate(curr_index-1), 'red should flash...') && strcmpi( LEDstate(curr_index-3), 'LEDTB On')
        issound(i) = 1;
        DelT(i)    = str2num( timestamps{ curr_index+1 } ) - str2num( timestamps{ curr_index-3+1 } );
        
    end;
    
    
end;
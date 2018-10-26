
function [isrepeat_same, isrepeat_different, norepeat] = parseRepeats(LEDstate)

s = strfind(LEDstate, 'autostart');
index = find(~cellfun(@isempty,s));
isrepeat = []; norepeat = []; success = [];

isrepeat_same          = zeros(1,length(index));
isrepeat_different     = zeros(1,length(index));
isrepeat               = zeros(1,length(index));
% isrepeat_same     = zeros(1,length(index));
norepeat          = zeros(1,length(index));

for i = 1:length(index)
    curr_index = index(i);
    if strcmpi( LEDstate{curr_index-3}, 'repeat cue same' )
        isrepeat_same(i) = 1;
    elseif strcmpi( LEDstate{curr_index-3}, 'repeat cue different' )
        isrepeat_different(i) = 1;
    elseif strcmpi( LEDstate{curr_index-3}, 'repeat cue' )
        isrepeat(i) = 1;
    else
        norepeat(i) = 1;
    end;
end;
%%

isrepeat_conflict = zeros(1, length(index));
isrepeat_same     = zeros(1, length(index));
indices_repeat    = index(find( isrepeat == 1));

for i = 1:length(isrepeat)
    if isrepeat(i) == 1
        
        curr_index =  index(i);
        
        if  ( strcmpi( LEDstate{curr_index-2}, 'LEDTB On' ) && strcmpi( LEDstate{curr_index-5}, 'LEDTG on') )...
                ||  ( strcmpi( LEDstate{curr_index-2}, 'LEDTB On' ) && strcmpi( LEDstate{curr_index-7}, 'LEDTG on') )
            isrepeat_different( i ) = 1;
            
            disp('Conflict');
        elseif  ( strcmpi( LEDstate{curr_index-2}, 'LEDTG On' ) && strcmpi( LEDstate{curr_index-5}, 'LEDTG on') )...
                ||  ( strcmpi( LEDstate{curr_index-2}, 'LEDTG On' ) && strcmpi( LEDstate{curr_index-7}, 'LEDTG on') )
            isrepeat_same( i ) = 1;
            
            disp('No Conflict');
        end;
    end;
end;


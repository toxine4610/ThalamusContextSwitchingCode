
function [isrepeat_same, isrepeat_different, norepeat] = findRepeats(LEDstate)


%%

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
    end
    
    if strcmpi( LEDstate(curr_index-1), 'red should flash...') && strcmpi( LEDstate(curr_index-4), 'repeat cue same')
        isrepeat_same(i) = 1;
    elseif strcmpi( LEDstate(curr_index-1), 'red should flash...') && strcmpi( LEDstate(curr_index-4), 'repeat cue different')
        isrepeat_different(i) = 1;
    end;
end;
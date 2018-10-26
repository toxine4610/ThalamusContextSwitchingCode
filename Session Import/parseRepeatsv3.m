
function [isrepeat_same, isrepeat_different, norepeat] = parseRepeatsv3(LEDstate)

s = strfind(LEDstate, 'autostart');
index = find(~cellfun(@isempty,s));
isrepeat = []; norepeat = []; success = [];

isrepeat_same          = zeros(1,length(index));
isrepeat_different     = zeros(1,length(index));
isrepeat               = zeros(1,length(index));
% isrepeat_same     = zeros(1,length(index));
norepeat          = zeros(1,length(index));

for i = 1:length(index)
    
    if i == 1
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'repeat cue different');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isrepeat_different(i-1) = 1;
        end
        clear ind
        
        
        s = strfind(foo, 'repeat cue same');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isrepeat_same(i) = 1;
        end
        clear ind
        
        
        s = strfind(foo, 'repeat cue');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isrepeat(i) = 1;
        end
        clear ind s
        
        s = strfind(foo, 'repeat cue');
        ind = find(~cellfun(@isempty,s));
        if isempty(ind);
            norepeat(i) = 1;
        end
        clear ind s
        
        
    elseif i >= 2
        currIndex_S = index(i-1);
        currIndex_E = index(i);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'repeat cue different');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isrepeat_different(i) = 1;
        end
        clear ind
        
        
        s = strfind(foo, 'repeat cue same');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isrepeat_same(i) = 1;
        end
        clear ind
        
        
        s = strfind(foo, 'repeat cue');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isrepeat(i) = 1;
        end
        clear ind s
        
        s = strfind(foo, 'repeat cue');
        ind = find(~cellfun(@isempty,s));
        if isempty(ind);
            norepeat(i) = 1;
        end
        clear ind s
        
        
        %     if strcmpi( LEDstate(currIndex_S:currIndex_E), 'repeat cue same' )
        %         isrepeat_same(i) = 1;
        %     elseif strcmpi( LEDstate(currIndex_S:currIndex_E), 'repeat cue different' )
        %         isrepeat_different(i) = 1;
        %     elseif strcmpi( LEDstate(currIndex_S:currIndex_E), 'repeat cue' )
        %         isrepeat(i) = 1;
        %     else
        %         norepeat(i) = 1;
        %     end;
        %
        
        
        %     if strcmpi( LEDstate{curr_index-3}, 'repeat cue same' )
        %         isrepeat_same(i) = 1;
        %     elseif strcmpi( LEDstate{curr_index-3}, 'repeat cue different' )
        %         isrepeat_different(i) = 1;
        %     elseif strcmpi( LEDstate{curr_index-3}, 'repeat cue' )
        %         isrepeat(i) = 1;
        %     else
        %         norepeat(i) = 1;
        %     end;
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
    

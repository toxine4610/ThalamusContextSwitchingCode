
function [isdim] = findDimTrials(LEDstate);

s = strfind(LEDstate, 'autostart');
index = find(~cellfun(@isempty,s));

isdim = zeros(1,length(index));

%------- Mix (search back...)
for i = 1:length(index)
    if i == 1
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'dim trial');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isdim(1) = 1;
        end
        clear ind
        
    elseif i >= 2
        currIndex_S = index(i-1);
        currIndex_E = index(i);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'dim trial');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isdim(i) = 1;
        end
        clear ind
        
    end;
end;
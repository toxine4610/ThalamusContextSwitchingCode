


function D = parseMouseYellow(fname)

[behavior] = textread(fname,'%s', -1,'delimiter','\t');

% %separating file into LED state and timestamps
% timestamps = behavior(1:2:end-1);
% 
% %convert to an aray of doubles
% LEDstate = behavior(2:2:end);


ct = 0;
ct2 = 0;

for i = 1:size(behavior,1)
    foo = sprintf( behavior{i} );
    v   = str2num(foo);
    if ~isempty(v)
        ct= ct+1;
        timestamps{ct} =  behavior{i};
    elseif isempty(v)
        ct2 = ct2 + 1;
        LEDstate{ct2} = behavior{i};
    end;
end;

LEDstate = LEDstate';
timestamps = timestamps';

%% sanity check
err = [];
for i = 1:length(LEDstate)
    foo = LEDstate{i};
    if ~isempty(str2num(foo))
        err = [err, i];
    end;
end;
if ~isempty(err)
    
    disp( LEDstate(err(1)) );
    error('Missing time stamp!!!');
else
    disp('No Errors');
end;


%%

s = strfind(LEDstate, 'autostart');
index = find(~cellfun(@isempty,s));

offsetcorrect = [];

if ~isempty(offsetcorrect)
    index = index( offsetcorrect:end);
end;

isdim      = zeros(1, length(index));
islaser    = zeros(1, length(index));

%%

%------- Dim Trials (search back...)
for i = 1:length(index)
    if i == 1
        
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'dim trial');
        ind = find(~cellfun(@isempty,s));
        
        if ~isempty(ind);
            isdim(i) = 1;
            x = sprintf(foo{ind(end)});
        end;
        
        clear ind ind2
        
    elseif i >= 2
        currIndex_S = index(i-1);
        currIndex_E = index(i);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'dim trial');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isdim(i) = 1;
            x = sprintf(foo{ind(end)});
        end;
        clear ind ind2
        
    end;
end;


%------- Laser Trials (search back...)
for i = 1:length(index)
    if i == 1
        
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'yellow laser on');
        ind = find(~cellfun(@isempty,s));
        
        if ~isempty(ind);
            islaser(i) = 1;
            x = sprintf(foo{ind(end)});
        end;
        
        clear ind ind2
        
    elseif i >= 2
        currIndex_S = index(i-1);
        currIndex_E = index(i);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'yellow laser on');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            islaser(i) = 1;
            x = sprintf(foo{ind(end)});
        end;
        clear ind ind2
        
    end;
end;

%%
indexDim     = find(isdim == 1);
indexLaser   = find(islaser == 1);

D  = zeros( length(index), 2);

D(indexDim,  1 ) = 1;
D(indexLaser,2) = 1;
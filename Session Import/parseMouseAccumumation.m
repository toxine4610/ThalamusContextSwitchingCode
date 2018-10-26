
function Out = parseMouseAccumumation(fname);

[behavior] = textread(fname,'%s', -1,'delimiter','\t');

% 
% %separating file into LED state and timestamps
% LEDstate = behavior(1:2:end-1);
% 
% %convert to an aray of doubles
% timestamps = behavior(2:2:end);

[LEDstate, timestamps] = correctMissingStamps(behavior);

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
% index = index(2:end);

isRule     = zeros(1,length(index));
iscorr     = zeros(1,length(index));
iswrong    = zeros(1,length(index));
ismix      = zeros(1,length(index));
mixvalue    = zeros(1, length(index));
CorrLocation    = zeros(1, length(index));
RespLocation    = zeros(1, length(index));
ruleType    = zeros(1, length(index));
isaborted  = zeros(1, length(index));
isdim      = zeros(1, length(index));

%% TimeStamps
Timestamp = [];

gapfind = 0;
gapOff = 0;

for i = 1:length(timestamps)
    timestamp = cell2mat(timestamps(i));
    timestamp = int32(str2num(timestamp));
    timestamp = double(timestamp);
    if (i > 2 && gapfind < 1)
        if timestamp < Timestamp(i-1)
            gapOff = Timestamp(i-1);
            gapfind = gapfind+1;
        end
    end
    timestamp = timestamp+gapOff;
    Timestamp = [Timestamp;timestamp];
end

auto_ind = find(strcmp(LEDstate, 'autostart') );
initT = Timestamp(auto_ind) - Timestamp(auto_ind - 1);

for n = 2:length(auto_ind)
    ITIval = Timestamp(auto_ind(n)) - Timestamp(auto_ind(n - 1));
    ITI(n) = ITIval;
end

for i = 1:length(index)
    autostart_ts(i) = str2num( timestamps{ index(i) } );
end;

Out.ats    = autostart_ts;
Out.initT = initT;
Out.ITI   = ITI;


%%
s2 = strfind(LEDstate, 'stim left on');
index2 = find(~cellfun(@isempty,s2));
for i = 1:length(index2)
    tsL(i) = str2num( timestamps{ index2(i) } );
end;

clear s2 index2
s2 = strfind(LEDstate, 'stim right on');
index2 = find(~cellfun(@isempty,s2));
for i = 1:length(index2)
    tsR(i) = str2num( timestamps{ index2(i) } );
end;

Out.ts = sort( [tsL,tsR] );

%% Behavior


% --- convention to be used here - Audio = 2 Outer pokes = 2, 
%                                - Visual = 1 Inner pokes = 1
%                                - Right= +, Left = -

%------ Removed
for i = 1:length(index)
    
    if ismember( i, [1:length(index)-1] )
        currIndex_S = index(i);
        currIndex_E = index(i+1);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'stim left on');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            foo2 = foo( ind:end );
            s = strfind(foo2, 'reward delivered');
            ind = find(~cellfun(@isempty,s));
            if isempty(ind)
                isaborted(i) = 1;
            end;
        end
        clear ind
    end;
    
    if i == length(index)
        currIndex_S = index(i);
        currIndex_E = length(LEDstate);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'removed');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            foo2 = foo( ind:end );
            s = strfind(foo2, 'reward delivered');
            ind = find(~cellfun(@isempty,s));
            if isempty(ind)
                isaborted(i) = 1;
            end;
        end
        clear ind
    end;
end;

for i = 1:length(index)
    
    if ismember( i, [1:length(index)-1] )
        currIndex_S = index(i);
        currIndex_E = index(i+1);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'stim right on');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            foo2 = foo( ind:end );
            s = strfind(foo2, 'reward delivered');
            ind = find(~cellfun(@isempty,s));
            if isempty(ind)
                isaborted(i) = 1;
            end;
        end
        clear ind
    end;
    
    if i == length(index)
        currIndex_S = index(i);
        currIndex_E = length(LEDstate);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'removed');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            foo2 = foo( ind:end );
            s = strfind(foo2, 'reward delivered');
            ind = find(~cellfun(@isempty,s));
            if isempty(ind);
                isaborted(i) = 1;
            end;
        end
        clear ind
    end;
end;



%------- Correct
for i = 1:length(index)
    
    if ismember( i, 1:length(index)-1 );
        currIndex_S = index(i);
        currIndex_E = index(i+1);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'reward delivered');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind)
            iscorr(i) = 1;
            x = sprintf(foo{ind(end)});
            if ~isempty( strfind( x, 'reward delivered OR') );
                cchoice(i) = 2;
            elseif ~isempty( strfind( x, 'reward delivered OL') );
                cchoice(i) = -2;
            elseif ~isempty( strfind( x, 'reward delivered IR') );
                cchoice(i)  = 1;
            elseif ~isempty( strfind( x, 'reward delivered IL') );
                cchoice(i) = -1;
            end;
        end
        clear ind
    end;
    
    if i == length(index)
        currIndex_S = index(i);
        currIndex_E = length(LEDstate);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'reward delivered');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            iscorr(i) = 1;
            x = sprintf(foo{ind(end)});
            if ~isempty( strfind( x, 'reward delivered OR') ); % 
                cchoice(i) = 2;
            elseif ~isempty( strfind( x, 'reward delivered OL') );
                cchoice(i) = -2;
            elseif ~isempty( strfind( x, 'reward delivered IR') ); 
                cchoice(i)  = 1;
            elseif ~isempty( strfind( x, 'reward delivered IL') );
                cchoice(i) = -1;
            end;
        end
        clear ind
    end;
end;



%------- response location
for i = 1:length(index)
    
    if ismember( i, 1:length(index)-1 );
        currIndex_S = index(i);
        currIndex_E = index(i+1);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'poked');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind)
            isResp(i) = 1;
            x = sprintf(foo{ind(end)});
            if ~isempty( strfind( x, 'poked Outer right') );
                RespLocation(i) = 2;
            elseif ~isempty( strfind( x, 'poked Outer left') );
                RespLocation(i) = -2;
            end;
        end
        clear ind
    end;
    
    if i == length(index)
        currIndex_S = index(i);
        currIndex_E = length(LEDstate);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'poked');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind)
            isResp(i) = 1;
            x = sprintf(foo{ind(end)});
            if ~isempty( strfind( x, 'poked Outer right') );
                RespLocation(i) = 2;
            elseif ~isempty( strfind( x, 'poked Outer left') );
                RespLocation(i) = -2;
            end;
        end
        clear ind
    end;
end;


%------- correct location
for i = 1:length(index)
    
    if ismember( i, 1:length(index)-1 );
        currIndex_S = index(i);
        currIndex_E = index(i+1);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'stim left on');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind)
            isChoice(i) = 1;
            CorrLocation(i) = -2;
        end
        clear ind
    end;
    
    if i == length(index)
        currIndex_S = index(i);
        currIndex_E = length(LEDstate);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'stim left on');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind)
            isChoice(i) = 1;
            CorrLocation(i) = -2;
        end
        clear ind
    end;
end;

for i = 1:length(index)
    
    if ismember( i, 1:length(index)-1 );
        currIndex_S = index(i);
        currIndex_E = index(i+1);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'stim right on');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind)
            isChoice(i) = 1;
            CorrLocation(i) = 2;
        end
        clear ind
    end;
    
    if i == length(index)
        currIndex_S = index(i);
        currIndex_E = length(LEDstate);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'stim right on');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind)
            isChoice(i) = 1;
            CorrLocation(i) = 2;
        end
        clear ind
    end;
end;


%------- RULE
for i = 1:length(index)
    
    if ismember( i, 1:length(index)-1 );
        currIndex_S = index(i);
        currIndex_E = index(i+1);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'audio stim');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isRule(i) = 1;
            ruleType(i) = 2;
        end
        clear ind
    end;
    
    if i == length(index)
        currIndex_S = index(i);
        currIndex_E = length(LEDstate);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'audio stim');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isRule(i) = 1;
            ruleType(i) = 2;
        end
        clear ind
    end;
end;

for i = 1:length(index)
    
    if ismember( i, 1:length(index)-1 );
        currIndex_S = index(i);
        currIndex_E = index(i+1);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'visual stim');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isRule(i) = 1;
            ruleType(i) = 1;
        end
        clear ind
    end;
    
    if i == length(index)
        currIndex_S = index(i);
        currIndex_E = length(LEDstate);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'visual stim');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isRule(i) = 1;
            ruleType(i) = 1;
        end
        clear ind
    end;
end;

%------- Mix (search back...)
for i = 1:length(index)
    if i == 1
        
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'cue seq');
        ind = find(~cellfun(@isempty,s));
        
        if ~isempty(ind);
            ismix(i) = 1;
            x = sprintf(foo{ind(end)});
            mixvalue(i) = str2double(x(9:12));
            
        elseif isempty(ind)
            clear s
            s = strfind(foo, 'cue seq');
            ind2 = find(~cellfun(@isempty,s));
            if ~isempty(ind2);
                ismix(i) = 1;
                x = sprintf(foo{ind2(end)});
                mixvalue(i) = str2double(x(9:12));
            end
        end;
        
        clear ind ind2
        
    elseif i >= 2
        currIndex_S = index(i-1);
        currIndex_E = index(i);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'cue seq');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            ismix(i) = 1;
            x = sprintf(foo{ind(end)});
            mixvalue(i) = str2double(x(9:12));
            
        elseif isempty(ind)
            clear s
            s = strfind(foo, 'cue seq');
            ind2 = find(~cellfun(@isempty,s));
            if ~isempty(ind2);
                ismix(i) = 1;
                x = sprintf(foo{ind2(end)});
                mixvalue(i) = str2double(x(9:12));
            end
        end;
        clear ind ind2
        
    end;
end;


%------- Dim Trials (search back...)
for i = 1:length(index)
    if i == 1
        
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'cue HP/LP');
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
        s = strfind(foo, 'cue HP/LP');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isdim(i) = 1;
            x = sprintf(foo{ind(end)});
        end;
        clear ind ind2
        
    end;
end;


for i = 1:length(index)
    if i == 1
        
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'cue uv/green');
        ind = find(~cellfun(@isempty,s));
        
        if ~isempty(ind);
            isdim(i) = 2;
            x = sprintf(foo{ind(end)});
        end;
        
        clear ind ind2
        
    elseif i >= 2
        currIndex_S = index(i-1);
        currIndex_E = index(i);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'cue uv/green');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isdim(i) = 2;
            x = sprintf(foo{ind(end)});
        end;
        clear ind ind2
        
    end;
end;


for i = 1:length(index)
    if i == 1
        
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'cue HP/green');
        ind = find(~cellfun(@isempty,s));
        
        if ~isempty(ind);
            isdim(i) = 3;
            x = sprintf(foo{ind(end)});
        end;
        
        clear ind ind2
        
    elseif i >= 2
        currIndex_S = index(i-1);
        currIndex_E = index(i);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'cue HP/green');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isdim(i) = 3;
            x = sprintf(foo{ind(end)});
        end;
        clear ind ind2
        
    end;
end;


for i = 1:length(index)
    if i == 1
        
        currIndex_S = index(i);
        foo = LEDstate(1:currIndex_S);
        s = strfind(foo, 'cue LP/uv');
        ind = find(~cellfun(@isempty,s));
        
        if ~isempty(ind);
            isdim(i) = 4;
            x = sprintf(foo{ind(end)});
        end;
        
        clear ind ind2
        
    elseif i >= 2
        currIndex_S = index(i-1);
        currIndex_E = index(i);
        
        foo = LEDstate(currIndex_S:currIndex_E);
        s = strfind(foo, 'cue LP/uv');
        ind = find(~cellfun(@isempty,s));
        if ~isempty(ind);
            isdim(i) = 4;
            x = sprintf(foo{ind(end)});
        end;
        clear ind ind2
        
    end;
end;

%%

iswrong( find( CorrLocation ~= RespLocation )) = 1;
isaborted( iswrong==1 ) = 0;

indexCorrect = find(iscorr==1);
indexWrong   = find(iswrong==1);
indexMix     = find(ismix == 1);


% build D matrix ........................................................
% col 1 = correct/wrong (1/0)
% col 2 = rule type (Vis - 1, Aud - 2)
% col 3 = smoothed performance
% col 4 = is it a Mixed trial (1/0)
% col 5 = what is the mix value
% col 6 = Correct target location
% col 7 = Response location
% col 8 = Context (1-HPLP, 2 - UV/Green, 3 - HP/Green, 4 -LP/UV)

D(indexCorrect, 1) = 1;
D(indexWrong, 1)   = 0;
D(:, 2) = ruleType;


w = gausswin(7);
w = w/sum(w);
D(:,3) = conv(D(:,1),w,'same');


D(indexMix, 4 ) = 1;
D(indexMix, 5 ) = 100.*(mixvalue);

D(:,6) = CorrLocation;
D(:,7) = RespLocation;
D(:, 8) = isdim;


Dcorrected = D( find(isaborted==0),:);
Dcorrected(:,2) = conv(Dcorrected(:,1),w,'same');

tscorrected = Out.ts( find(isaborted==0) );
Out.tscorrected = tscorrected;
Out.initTcorrected = Out.initT( find(isaborted==0) );
Out.atscorrcted = Out.ats( find(isaborted==0) );

Out.indexRemoved = find(isaborted == 1);
Out.indexOK      = find(isaborted == 0);

Out.D = D;
Out.Dcorrected = Dcorrected;

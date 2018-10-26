

function [raster_data,raster_labels,raster_site_info] = Mixed_prepDateforDecodingToolbox(foo,session_num);

ID = { 'R1_sound', 'R2_sound', 'R1_nosound', 'R2_nosound' };

Identifiers =  {'BrownNoise', 'BlueNoise', 'GreenLED', 'BlueLED' };
Context      = {'AudCue', 'AudCue', 'LEDCue', 'LEDCue' };
Descriptions = {'A2V', 'A2A', 'A2V', 'A2A'};



%% build time vector

endPoints = [-0.2 1.2];

timeSize = (endPoints(2)-endPoints(1))*1000;
timeVector = endPoints(1):0.001:endPoints(2);


%%

clear raster_labels raster_site_info
raster_data = [];
stimulus_ID = {};
stimulus_position = {};
stimulus_context = {};
combined_ID_position = {};
combined = [];

for n = 1:4
    
    if ismember(n,1:2);
        cellModulation = 'Context1';
    elseif ismember(n,3:4);
        cellModulation = 'Context2';
    end;
    
    
    zeroMatrix = zeros(1,timeSize+1);
    
    eval(['currentDataSingleType = foo.SpikeTimes_' ID{n} ';'])
    
    for nn = 1:length(currentDataSingleType)
        currentDataSingleTypeCell = currentDataSingleType{nn};
        holderRow = histc(currentDataSingleTypeCell,timeVector)';
        if length(holderRow > 0)
            zeroMatrix(nn,:) = histc(currentDataSingleTypeCell,timeVector)';
        else
            zeroMatrix(nn,:) = 0;
        end
        currStimID{nn}       = Identifiers{n};
        currTrialContext{nn} = Context{n};
        currRuleMeaning{nn}  = [Descriptions{n}];
        currCombined{nn}     = [Descriptions{n} '_' Identifiers{n}];
        %         currStimContext{nn} = Context{n};
    end
    
    raster_data = [raster_data;zeroMatrix];
    %     stimulus_context = [stimulus_context; currStimContext];
    stimulus_ID = [stimulus_ID; currStimID'];
    stimulus_position = [stimulus_position; currTrialContext'];
    combined_ID_position = [combined_ID_position;   currRuleMeaning'];
    combined  = [combined; currCombined'];
    clear zeroMatrix currStimID currTrialContext currRuleMeaning currCombined;
    
end

raster_labels.stimulus_ID   = stimulus_ID;
raster_labels.trial_context = stimulus_position;
raster_labels.rule_meaning  = combined_ID_position;
raster_labels.cominedInfo   = combined;

raster_site_info.session_ID = session_num;
% raster_site_info.unit = num2str(foo.unitNbr);
raster_site_info.alignment_event_time = 1*1000;




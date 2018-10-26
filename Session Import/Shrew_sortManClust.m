

% close all
pre = 5;
post = 3;
clearTrue = 0;
lasShowVal = 1;
baseEnds = [-0.5 0];
filtSize = [5 2];
postTag = 'baseline';

includeSet = [1:size(num_seq,1)];

% pure correct
VisL_corr = find( D(:,1) == 1 & D(:,5) == 1 );
VisR_corr = find( D(:,1) == 1 &  D(:,5) == -1 );
Vis_corr  = sort( [VisL_corr; VisR_corr] );

% pure correct
AudL_corr = find( D(:,1) == 1  & D(:,5) == 2 );
AudR_corr = find( D(:,1) == 1  & D(:,5) == -2 );
Aud_corr  = sort( [AudL_corr; AudR_corr] );


% pure incorrect
VisL_incorr = find( D(:,1) == 0  & D(:,5) == 1 );
VisR_incorr = find( D(:,1) == 0  & D(:,5) == -1 );
Vis_incorr  = sort( [VisL_incorr; VisR_incorr] );

% pure incorrect
AudL_incorr = find( D(:,1) == 0  & D(:,5) == 2 );
AudR_incorr = find( D(:,1) == 0  & D(:,5) == -2 );
Aud_incorr  = sort( [AudL_incorr; AudR_incorr] );


for nn = 1:size(num_seq,1)
    
    for m = 1:4
        
        if m == 1
            f = Vis_corr;
            stringVal = 'Vis_corr';
            stringVals{m} = stringVal;
        elseif m == 2
            f = Aud_corr;
            stringVal = 'Aud_corr';
            stringVals{m} = stringVal;
            
            
        elseif m == 3
            f = Vis_incorr;
            stringVal = 'Vis_incorr';
            stringVals{m} = stringVal;
        elseif m == 4
            f = Aud_incorr;
            stringVal = 'Aud_incorr';
            stringVals{m} = stringVal;
            
        end
        eval(['cl' num2str(nn) '_align = align(cl' num2str(nn) ',Linit,f,pre,post);']);
        eval(['clcurr_align = cl' num2str(nn) '_align;'])
        eval(['clcurr_align' num2str(m) ' = cl' num2str(nn) '_align;']) % save the raster spike timing CG
        
        timestamps = clcurr_align.timestamp;
        trial_id = clcurr_align.trial_id';
        trials = unique(trial_id);
        
        clear clCell_align
        clCell_align = []
        qAdd = 1;
        
        for q = 1:length(f);
            if ~isempty(intersect(q,trials))
                clCell_align{q} = timestamps(find(trial_id == q));
                %                 clCell_align{q} = clCell_align{q}(intersect(find(clCell_align{q} > winS),find(clCell_align{q} < winEnd)));
            else
                clCell_align{q} = [];
            end
            fine = f;
        end
        fine = 1:length(clCell_align);
        fine = fine';
        eval(['cl' num2str(nn) '_cellData' postTag '.' stringVal ' = clCell_align;'])
    end

end;
function eyeposF = upsampleEyepos(eyepos, Fs)
% up-sample and spline interpolate eye position to sampling rate
% eyeposF = upsampleEyepos(eyepos, Fs)
% eyepos is cell-array of [time x y] eye position vectors

if nargin < 2
    Fs = 1e3;
end

bs = 1/Fs;

nTrials = numel(eyepos);

eyeposF = cell(nTrials, 1);

for kTrial = 1:nTrials
    
    start = min(eyepos{kTrial}(:,1));
    stop  = max(eyepos{kTrial}(:,1));
    
    time_current = eyepos{kTrial}(:,1);
    
    time_new  = (start:bs:stop)';
    nsamples  = numel(time_new);
    x_new     = nan(nsamples, 1);
    y_new     = nan(nsamples, 1);
    
    
    dt = abs(bsxfun(@minus, time_new',time_current));
    
    [~, idx] = min(dt, [], 2);
    x_new(idx) = eyepos{kTrial}(:,2);
    y_new(idx) = eyepos{kTrial}(:,3);
    
    dx = mean([max(diff(find(diff(isnan(x_new))==1))) mean(diff(find(diff(isnan(x_new))==1)))]);
    dy = mean([max(diff(find(diff(isnan(x_new))==1))) mean(diff(find(diff(isnan(x_new))==1)))]);
    
    
    if isnan(dx) || isnan(dy)
        continue
    end
    
    dx = max(dx, 14);
    dy = max(dy, 14);
    % upsample and smooth
    x_filt = filtfilt(ones(ceil(dx),1)/ceil(dx), 1, repnan(x_new, 'spline'));
    y_filt = filtfilt(ones(ceil(dy),1)/ceil(dy), 1, repnan(y_new, 'spline'));
    
    eyeposF{kTrial} = [time_new(:) x_filt(:) y_filt(:)];
    
end
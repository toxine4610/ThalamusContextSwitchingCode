function [w,widxs, widx]=binContinuous(an,ev,win)
% bin sampled signal
% [w,widxs, widx]=binContinuous(an,ev,win)
% an is a column vector!!
if size(an,1)<size(an,2)
    an=an';
end
widx = win(1):win(2);

% ev(isnan(ev)) = []; % remove impossible spikes
% ev(ev > size(an,1)) = [];

widxs = bsxfun(@plus, ev(:), widx);
widxs(widxs(:,1) <= 0,:) = nan; % ignore very early spikes
widxs(widxs>size(an,1))=nan;
% widxs(widxs(:,end) > size(an,1),:) = nan; % ignore very late spikes

w=nan(size(widxs));
w(~isnan(widxs))=an(widxs(~isnan(widxs)));
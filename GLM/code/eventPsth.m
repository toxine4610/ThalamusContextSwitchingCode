function [m,s,bc,v, tspcntOut] = eventPsth(sptimes, ev, win, bs, skern)
% Make Peri-Stimulus-Time-Histogram (PSTH) aligned to stimulus events
% [m, s, bc, v, tspcnt] = eventPsth(sptimes, ev, win, bs, skern)
%   inputs
%       sptimes [n x 1] - vector of spike times
%       ev      [m x 1] - vector of event times
%       win     [1 x 2] - window around events to count spikes in
%       bs      [1 x 1] - scalar bin size
%       skern   [1 x k] - smoothing kernel
%       works in any units as long as all variable have the same units.

% 20140415 jly wrote it


if nargin < 5
    skern = [];
end

ev=ev(:);

be = win(1):bs:win(2);
bc = be(1:end-1)+bs/2;

% normalize smoothing kernel
if ~isempty(skern)
    k = numel(skern);
    wins=[(win(1)-bs*k) (win(2)+bs*k)];
    be = (win(1)-bs*k):bs:(win(2)+bs*k);
    bc = be(1:end-1)+bs/2;
    goodindex = bc>win(1) & bc <win(2);
    skern = skern/sum(skern);
else
    wins=win;
end

[tspcnt,~,validEvents]=binSpTimes(sptimes, ev, wins, bs);

% if smoothing needs to be done
if ~isempty(skern)
    %     tspcnt = filter(skern, 1, tspcnt, [], 2);
    %     nz = filter(skern, 1, ones(size(tspcnt)), [], 2);    
    tspcnt(validEvents,:) = filtfilt(skern, 1, tspcnt(validEvents,:)')';
    nz = filtfilt(skern, 1, ones(size(tspcnt))')';
    tspcnt(validEvents,:) = tspcnt(validEvents,:)./nz(validEvents,:);
    tspcnt = tspcnt(:, goodindex);
    bc = bc(goodindex);
end

% get summary statistics
m = nanmean(tspcnt)/bs;
v = nanvar(tspcnt)/bs;
s = nanstd(tspcnt)/sqrt(sum(validEvents))/bs;

tspcntOut = tspcnt;
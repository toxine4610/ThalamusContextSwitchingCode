function [rbar, psthTime, rse]=computeChoPSTH(y, g, trial, aligningField, win)
% compute the coherence dpendent psth using spike predictions and mstruct
% [rbar, psthTime]=computeChoPSTH(y, g, trial)
if exist('win', 'var')
    win=win/g.binSize;
else
    win=[-10 150];
end

if nargin<4 || isempty(aligningField)
    aligningField='motionon';
end
   

motionOn=find(g.getBinnedSpikeTrain(trial, aligningField, 1:numel(trial)));
cho=arrayfun(@(x) x.choice, trial);


psthTime=(win(1):win(2))*g.binSize;
[~,~,~,rtrue]=eventTriggeredAverage(y, motionOn, win); %#ok<FNDSB>

% dpcell=dprime(nanmean(rtrue(:,psthTime>0 & psthTime < 1e3),2),cho);
% if dpcell<0
%     cho=~cho;
% end
            
rbar=[nanmean(rtrue(cho,:))' nanmean(rtrue(~cho,:))'];
rse=[(nanstd(rtrue(cho,:))')/sqrt(sum(cho)) (nanstd(rtrue(~cho,:))')/sqrt(sum(~cho))];

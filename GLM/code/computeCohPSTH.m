function [rbar, psthTime, nTrials, q, dpsign]=computeCohPSTH(y, g, trial, win, plotIt, choFlip, dpcell)
% compute the coherence dpendent psth using spike predictions and mstruct
% [rbar, psthTime, nTrialsPerCond, Conditions]=computeCohPSTH(y, g, trial, win, plotIt)
if ~exist('win', 'var') || isempty(win)
    win=[-10 120];
end
    
if ~exist('choFlip', 'var')
    choFlip=true; % flip to reflect the prefered direction
end

if ~exist('plotIt', 'var')
    plotIt=false;
end

if isfield(trial, 'trial')
    trial=trial.trial;
end

motionOn=find(g.getBinnedSpikeTrain(trial, 'motionon', 1:numel(trial)));
% pulses=cell2mat(arrayfun(@(x) x.pulses, trial(:), 'UniformOutput', false));
% pulses=pulses/std(pulses(:));
nTrials=numel(trial);
coh=nan(nTrials,1);
for kTrial=1:nTrials
coh(kTrial)=sum(trial(kTrial).pulses);
end
% coh=mean(pulses,2);
% coh=zscore(coh);
% coh=zscore(pulses(:,1));

cho=arrayfun(@(x) x.choice, trial);
% if mean((coh>0)==cho)<.5
%     coh=-coh;
% end
% 
% p=pmfTools(abs(coh), (coh>0)==cho);
% % figure, 
% % p.plot
% p.fit
% p.getThreshold([.6 .8])


psthTime=(win(1):win(2))*g.binSize;
nPsthBins=numel(psthTime);
[~,~,~,rtrue]=eventTriggeredAverage(y, motionOn, win); %#ok<FNDSB>
% [rtrue, bcenters, vidx]=binSpTimes(y
if ~exist('dpcell', 'var')
    ix=sum(isnan(rtrue(:,psthTime>100 & psthTime < 1200)),2)==0;
    rtr=rtrue(ix,psthTime>0 & psthTime < 1400);
    dpcell=dprime(nanmean(rtr,2),coh(ix)>0);
end

if choFlip
    dpsign=dpcell;
else
    dpsign=1;
end

if dpcell<0 && choFlip
    cho=~cho;
end
            
abscoh=abs(coh);
abscoh(abscoh==0)=min(abscoh(abscoh~=0));
scoh=sign(cho-.5).*abscoh; % aligned to choice now


%% use coherence sorting
% for right choices

binId=nan(nTrials,1);

rix=find(scoh>0);
[s, id]=sort(scoh(rix));
n=numel(s);
step=ceil(n/4);

for k=1:4
    ii=(1+(k-1)*step):(k*step);
    if k==4
        ii=(3*step):n;
    end
    binId(rix(id(ii)))=k;
end

lix=find(scoh<0);
[s, id]=sort(scoh(lix));
n=numel(s);
step=ceil(n/4);
for k=1:4
     ii=(1+(k-1)*step):(k*step);
    if k==4
        ii=(3*step):n;
    end
    binId(lix(id(ii)))=k+4;
end

nCohs=numel(unique(binId));
rbar=nan(nPsthBins, nCohs);
nTrials=nan(1,nCohs);
q=nan(1,nCohs);
for k=1:nCohs
    rbar(:,k)=nanmean(rtrue(binId==k,:));
    nTrials(k)=sum(binId==k);
    q(k)=mean(scoh(binId==k));
end

rbar(:,1:4)=fliplr(rbar(:,1:4));
rbar(:,5:end)=fliplr(rbar(:,5:end));
q=q([4 3 2 1 8 7 6 5]);
nTrials=nTrials([4 3 2 1 8 7 6 5]);
rbar=fliplr(rbar);
q=fliplr(q);
nTrials=fliplr(nTrials);
%%

% %% OLD WAY
% q=quantile(scoh(scoh>0), [.25 .5 .75 1]);%[.33 .66 1]
% mq=quantile(scoh(scoh<0), [0 .25 .5 .75]);
% % q=[p.threshValue inf];
% binEdges=[mq 0 q];
% % binEdges = sort([-q 0 q]);
% 
% % binEdges = q;
% nCohs=numel(binEdges)-1;
% if isfield(trial, 'trialCnt')
%     goodix=[trial(:).trialCnt]<10;
% else
%     goodix=true(nTrials,1);
% end
% 
% binId=cell2mat(arrayfun(@(x,y) goodix(:) & scoh>=x & scoh<=y, binEdges(1:end-1), binEdges(2:end), 'UniformOutput', false));
% 
% 
%             
% rbar=nan(nPsthBins, nCohs);
%             
% for k=1:nCohs
%         rbar(:,k)=nanmean(rtrue(binId(:,k),:));
% end
% 
% nTrials=sum(binId);
%     

if nargout<1 || plotIt
    set(gcf, 'defaultAxesColorOrder', cbrewer('jake', 'rdbu', nCohs))
    plot(psthTime, rbar)
    xlim([-100 1200])
end
%%
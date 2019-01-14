function [rbar,ptaTime, fl, dpcell]=computePTA(y, g, trial, choFlag, win, nk, dpcell, idx)
% compute the whitened PTA using spike train and mstruct
% [rbar,ptaTime, fl]=computePTA(y, g, trial, choFlag, win, nk)
if ~exist('win', 'var')
    win=ceil([-100 1600]./g.binSize);
end

if ~exist('choFlag', 'var')
    choFlag=true;
end

if ~exist('nk' ,'var')
    nk=60;
end

if isfield(trial, 'trial')
    trial=trial.trial;
end

if ~exist('idx', 'var')
    idx = true(numel(trial),1);
end
    
if ~islogical(idx)
    inds = idx;
    idx = false(numel(trial),1);
    idx(inds) = true;
end

    
%     b=neuroGLM('binSize', 10);
%     nk=nk/10;
%     win=win./10;
%     dur=numel(b.getBinnedSpikeTrain(trial, 'motionon', 1:numel(trial)));
%     dur2=numel(g.getBinnedSpikeTrain(trial, 'motionon', 1:numel(trial)));
%     if issparse(y)
%         sptimes=find(y); % only one spike per bin
%         bsptimes=b.binfun(sptimes); %#ok<FNDSB>
%         y=full(sparse(bsptimes, 1, 1, dur, 1));
%     else
%         y0=y(1:10:end);
%     end
% else
%     b=neuroGLM('binSize', g.binSize);
%     
% end
% motionOn=find(g.getBinnedSpikeTrain(trial, 'motionon', g.dm.trialIndices));
motionOn=g.getBinnedSpikeTrain(trial, 'motionon', 1:numel(trial));

ptimes=g.binfun(-g.binSize*win(1)+trial(1).pulseon-trial(1).pulseon(1));
ptaTime=(win(1):win(2))*g.binSize;

if g.binSize==1 % too small. make it 10ms bins
   win=win/10;
   nk=nk/10;
   motionOn=smooth(full(motionOn), 11);
   motionOn=motionOn(1:10:end);
   y=smooth(y,10)*10;
   y=y(1:10:end);
   ptaTime=ptaTime(1:10:end);
   ptimes=ptimes/10;
end

motionOn=find(motionOn);
motionOn(diff(motionOn)==1)=[];


badix=(numel(y)-motionOn)<win(2);
badix = badix | ~idx;

rbar = nan(nk, 7);
fl = nan;
if sum(~badix) < 25 
    dpcell = nan;
    return
end

cho=arrayfun(@(x) x.choice, trial(~badix));


rtrue=binContinuous(y, motionOn(~badix), win);


pulses=cell2mat(arrayfun(@(x) x.pulses, trial(~badix), 'UniformOutput', false)); %*trialParam.nGabors;
[pulses,mu, sigma]=zscore(pulses);
% mu=0; sigma=1;

%%
if ~exist('dpcell', 'var')
    dpcell=dprime(nanmean(rtrue(:,ptaTime>0 & ptaTime < 1e3),2),cho);
end

if dpcell<0
    pulses=-pulses;
end

r=sum(rtrue(:,ptaTime>0 & ptaTime < 1e3),2);

if g.binSize==1
    ptaTime=bsxfun(@plus, ptaTime((1:nk))', 10*ptimes(:)');
else
    ptaTime=bsxfun(@plus, ptaTime((1:nk))', g.binSize*ptimes(:)');
end
rc=rank(cov(pulses));
% fprintf('rank of pulses is %d\n', rc)
% if mean(r) <5 || rc < 7
%     rbar=nan(size(ptaTime));
%     return
% end
fl=[mean(r) rc];
%
switch choFlag
    case 1
        rcho=rtrue;
        rcho(cho,:)=ones(sum(cho),1)*mean(rcho(cho,:));
        rcho(~cho,:)=ones(sum(~cho),1)*mean(rcho(~cho,:));
        %     rbar=pulseSTA(pulses, rcho, ptimes, nk, 1);
        rbar=pulseSTA(pulses, rtrue-rcho, ptimes, nk, 1);
    case 2
        rcho=rtrue;
        rcho(cho,:)=ones(sum(cho),1)*mean(rcho(cho,:));
        rcho(~cho,:)=ones(sum(~cho),1)*mean(rcho(~cho,:));
        rbar=pulseSTA(pulses, rcho, ptimes, nk, 1);
    otherwise
        rbar=pulseSTA(pulses, bsxfun(@minus, rtrue, mean(rtrue)), ptimes, nk, 1);
end
rbar=bsxfun(@rdivide, rbar, sigma);
rbar=bsxfun(@minus, rbar, mu);


function [spcnt, bcenters,validEvents]=binSpTimes(sptimes, ev, win, binSize)
% bin spike times around event
% [spcnt, bcenters,validEvents]=binSpTimes(sptimes, ev, win, binSize)


if nargin < 4
    dsptimes=diff(sptimes);
    if nargin < 3
        win = [-mean(dsptimes) 10*mean(dsptimes)];
        if nargin < 2
            help binSpCounts
            return
        end
    end
    binSize = min(dsptimes);
end
sptimes=sptimes(:);

be = win(1):binSize:win(2);
bcenters = be(1:end-1)+binSize/2;

ev=ev(:);
binfun = @(t) (t == 0) + ceil(t/binSize);

validEvents = ~isnan(ev);
nTrials=numel(ev);

nEvents = sum(validEvents);

sbn = [];
str = [];
% assert(nEvents<2e3, 'too many events for this to run fast!')
for kEvent=1:nTrials
    if ~validEvents(kEvent)
        continue
    end
    spo = sptimes(sptimes > ev(kEvent) + be(1) & sptimes < ev(kEvent) + be(end))- ev(kEvent);
    if ~isempty(spo)
    sbn = [sbn; binfun(spo- be(1))]; %#ok<AGROW>
    str = [str; ones(numel(spo),1)*kEvent]; %#ok<AGROW>
    end
end

spcnt = full(sparse(str, sbn, 1, nTrials, binfun(be(end)-be(1))));
spcnt(~validEvents,:)=nan;
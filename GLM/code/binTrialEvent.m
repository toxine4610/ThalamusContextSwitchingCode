function cnt=binTrialEvent(trial, event, aligningField, win, binfun)
% bin a trialevent around an aligning field
% cnt=binTrialEvent(trial, event, aligningField, win, binfun)

if numel(trial(1).(event))==trial(1).duration
    isContinuous=true;
else
    isContinuous=false;
end

assert(isfield(trial(1), aligningField), 'aligningField must be a field of trial')
assert(isfield(trial(1), event), 'event must be a field of trial')
assert(isa(binfun, 'function_handle'), 'binfun must be a function handle to bin events')
nTrials=numel(trial);
nBins=binfun(diff(win));

if isContinuous
    % find bin size
    binSize=1;
    while binfun(binSize)==1
        binSize=binSize+1;
    end
    binSize=binSize-1;
    [~,binMat]=basisFactory.binDataVector(zeros(diff(win),1), binSize);
    nBins=binfun(diff(win));
    cnt=zeros(nTrials,nBins);
    for kTrial=1:nTrials
        t=trial(kTrial).(aligningField);
        ix=(win(1):win(2)-1)+t;
        validIx=ix>0&ix<=trial(kTrial).duration;
        trx=trial(kTrial).(event)(ix(validIx));
        cnt(kTrial,:)=trx(:)'*binMat(1:sum(validIx),:);
    end

else
    st=[];
    sts=[];
    for kTrial=1:nTrials
        x=trial(kTrial);
        tst=binfun(x.(event)(x.(event) > (win(1)+x.(aligningField)) & (x.(event) < (win(2)+x.(aligningField))))-x.(aligningField));
        sts=[sts; tst(:)]; %#ok<*AGROW>
        st=[st; ones(numel(tst),1)*kTrial];
    end
    stimes=sts-binfun(win(1));
    st(stimes<=0)=[];
    stimes(stimes<=0)=[];
    cnt=full(sparse(st, stimes, 1, nTrials, nBins));
end
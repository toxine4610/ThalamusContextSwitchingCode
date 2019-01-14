function PTA = pulseSTASingle(pulses, rtrial, ptimes, nkt)
% PTA = pulseSTA(pulses, rtrial, ptimes, nkt)

nb = size(rtrial,2);


nTrial = size(rtrial,1);


%% New crazy hack for speed
if nkt>size(rtrial,2)
    PTA=nan(nkt, 1);
    return
end

%%
tstarts=(0:nTrial-1)*nb;
t=bsxfun(@plus, ptimes, tstarts);
t=t(:);
p=reshape(pulses', [], 1);
n=numel(t);
X=sparse(t, ones(n,1), p, nb*nTrial, 1);
Xs=makeStimRows(X, nkt);
% y=reshape(rtrial', [], 1);
% % 
% % X0=ones(nb*nTrial,1);
% % Xs=[X0 Xs];
% % 
% % y=reshape(rtrial', [], 1);
% % y=smooth(y,10);
% % PTA=(Xs'*Xs)\(Xs'*y);
% % PTA(1)=[];
% % PTA=flipud(PTA);
% % plot(PTA)
% % 
% % %%
% % y=reshape(rtrial', [], 1);
% % 
% % % y=smooth(y,10);
% % y=y-mean(y);
% % 
% % Xs=makeStimRows(y, nkt);
% % X=sparse(t, repmat((1:7)', nTrial,1), p, nb*nTrial, 7);
% % 
% % Xs=Xs(nkt:end,:);
% % Xs=[Xs ones(size(Xs,1),1)];
% % X=sum(X(1:end-(nkt-1),:),2);
% % % X=X(1:end-(nkt-1),:);
% % % (Xs'*Xs + 1e7*eye(size(Xs,2)))\
% % PTA=(Xs'*X)/sum(X); %(Xs'*Xs)\
% % plot(PTA)
% % 
% % %%
% % % nPulses = size(pulses,2);
% % %     xi=1:nPulses:n;
% % %     X0=sparse(t(xi), ones(numel(xi),1), ones(numel(xi),1), nb*nTrial,1);
% % %     X0=makeStimRows(X0, nkt*2);


X0=ones(nb*nTrial,1);
Xs=[X0 Xs];

y=reshape(rtrial', [], 1);
PTA=(Xs'*Xs)\(Xs'*y);
PTA(1)=[];
PTA=flipud(PTA);
%     %%
%     figure(1); clf
%
%     plot(r);
%
%
%     %%
%
%     nY = numel(rtrial);
%     if flag
%
%         PTA = (reshape(XX\XY, [], nPulses+nc));
%         PTA(:,1:nc) = []; % remove bias
% %         PTA = reshape(XX(nkt+1:end, nkt+1:end)\XY(nkt+1:end), [], nPulses);
%     else
%         xdat.xx = XX;
%         xdat.xy = XY;
%         xdat.yy = YY;
%         xdat.ny = nY;
%
%         colStartInds = (0:nkt:((nPulses+nc)*nkt))+1;
%         khat = autoCorrRidgeRegress(xdat, colStartInds);
%         PTA  = (reshape(khat, [], nPulses+nc));
%         PTA(:,1:nc) = []; % remove bias
%     end
%
%
% if ~isempty(choices)
%     PTA(:,1:2) = [];
% end

% if sum(PTA(:,1))<0
%     PTA = -PTA;
% end
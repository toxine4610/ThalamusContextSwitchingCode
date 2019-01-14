function PTA = pulseSTA(pulses, rtrial, ptimes, nkt, flag, choices)
% PTA = pulseSTA(pulses, rtrial, ptimes, nkt)

if ~exist('flag', 'var')
    flag = false;
end

if ~exist('choices', 'var')
    choices = [];
end

nb = size(rtrial,2);
nPulses = size(pulses,2);

nTrial = size(rtrial,1);


%% New crazy hack for speed
if nkt>size(rtrial,2)
    PTA=nan(nkt, 7);
    return
end
try
    if isempty(choices)
        nc=1;
        filtk=diag(ones(nkt,1));
        kFilt=[zeros(ptimes(1)-1, nkt); filtk; zeros(nb-nkt-ptimes(1)+1, nkt)];
        Xs=kFilt;
        XY=sum(rtrial*kFilt)';
        for kPulse=1:nPulses
            kFilt=[zeros(ptimes(kPulse)-1, nkt); filtk; zeros(nb-nkt-ptimes(kPulse)+1, nkt)];
            Xs=[Xs kFilt];
            XY=[XY; sum(bsxfun(@times, (rtrial*kFilt), pulses(:,kPulse)))'];
        end
        
        
        YY=0;
        XX=0;
        for kTrial=1:nTrial
            XX=XX+[1 pulses(kTrial,:)]'*[1 pulses(kTrial,:)];
            YY = YY + rtrial(kTrial,:)*rtrial(kTrial,:)';
        end
        XX=(Xs'*Xs).*imresize(XX, nkt, 'nearest');
        
    else
        
        %% old slow (readable) way
        XX = 0;
        XY = 0;
        YY = 0;
        X = zeros(nb, nPulses);
        for kPulse = 1:nPulses
            X(ptimes(kPulse), kPulse) = 1;
        end
        X = [X(:,1) X]; % column for motion onset
        nc = 1;
        if ~isempty(choices)
            X = [X(:,1) X(:,1) X];
            nc = nc + 2;
        end
        
        Xs = makeStimRows(X, nkt);
        
        % xxt=Xs'*Xs; % precompute covariance of timing
        
        for kTrial = 1:nTrial
            if ~isempty(choices)
                Xb = bsxfun(@times, Xs, [ones(1,nkt) choices(kTrial)*ones(1,nkt) ~choices(kTrial)*ones(1,nkt) reshape(repmat(pulses(kTrial,:), nkt, 1), [], 1)']);
            else
                
                Xb = bsxfun(@times, Xs, [ones(1,nkt) reshape(repmat(pulses(kTrial,:), nkt, 1), [], 1)']);
                
            end
            
            % slightly faster version
            % 	xxp=imresize([1 pulses(kTrial,:)]'*[1 pulses(kTrial,:)], nkt, 'nearest');
            %     XX = XX + xxp.*xxt;
            % Slow version
            XX = XX + Xb'*Xb;
            
            XY = XY + Xb'*rtrial(kTrial,:)';
            YY = YY + rtrial(kTrial,:)*rtrial(kTrial,:)';
            
        end
    end
    
    %%
    
    nY = numel(rtrial);
    if flag
        
        PTA = (reshape(XX\XY, [], nPulses+nc));
        PTA(:,1:nc) = []; % remove bias
%         PTA = reshape(XX(nkt+1:end, nkt+1:end)\XY(nkt+1:end), [], nPulses);
    else
        xdat.xx = XX;
        xdat.xy = XY;
        xdat.yy = YY;
        xdat.ny = nY;
        
        colStartInds = (0:nkt:((nPulses+nc)*nkt))+1;
        khat = autoCorrRidgeRegress(xdat, colStartInds);
        PTA  = (reshape(khat, [], nPulses+nc));
        PTA(:,1:nc) = []; % remove bias
    end
    
    
catch me
    PTA=nan(nkt, 7);
end
% if ~isempty(choices)
%     PTA(:,1:2) = [];
% end

% if sum(PTA(:,1))<0
%     PTA = -PTA;
% end
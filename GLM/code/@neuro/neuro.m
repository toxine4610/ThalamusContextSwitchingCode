classdef neuro < handle
    properties
        exname
        nreader
        directory
        filename
        stim
    end
    
    methods (Access=public)
        function N = neuro(filename) % constructor
            if exist(filename, 'file')
                [pathto, fname, ext] = fileparts(filename);
                N.directory = pathto;
                N.filename  = [fname ext];
                N.nreader   = reader(fullfile(N.directory,N.filename));
                N.exname    = N.nreader.exname;
            else
               error('That file does not exist')
            end
            N.stim = rcgreader(N.exname, [N.directory '/../']);
        end
        
        %--- get neuron struct
        function S = getStruct(N)
            if exist(fullfile(N.directory, N.filename), 'file')
                S = load(fullfile(N.directory, N.filename));
            end
        end
        
        %--- which monkey?
        function m = monkey(N)
            switch N.exname(1)
                case 'p'
                    m = 1;
                case 'n'
                    m = 2;
            end
        end
        
        
        %--- Distance from stimulus tuning
        function p = deltaStim(N)
            p = diffCirc(N.prefDir,N.stim.theta);
        end
        
        
        %--- Is an MT neuron
        function b = isMT(N)
            b = strcmp(N.nreader.brainArea, 'MT');
        end
        
        %--- Is an LIP neuron
        function b = isLIP(N)
            b = strcmp(N.nreader.brainArea, 'LIP');
        end
        
        %--- nTrials
        function n = nTrials(N)
            n = numel(N.nreader.trialIndex);
        end
        
        
        function d=plotTaskGeometry(N)
            d=gca;
            cla
            pos=N.stim.gaborXY';
            nGabors=size(pos,2);
            sigma=.25; %N.stim.sf;
            % set up stimulus space
            mnx=-25; %min(pos(1,:))-min(abs(diff(pos(1,:))));
            mny=-15; %min(pos(2,:))-min(abs(diff(pos(1,:))));
            mxx=25; %max(pos(1,:))+min(abs(diff(pos(1,:))));
            mxy=15; %max(pos(2,:))+min(abs(diff(pos(1,:))));
            
            xax=mnx:(1/N.stim.ppd):mxx;
            yax=mny:(1/N.stim.ppd):mxy;
            
            [xx,yy]=meshgrid(xax,yax);
            Cframe=zeros(size(xx));
            
            % build spatial ramp for orientation of sinewave carrier
            ramp=cosd(N.stim.theta)*xx + sind(N.stim.theta)*yy;
            
            for kGabor=1:nGabors
                g=exp(-((xx-pos(1,kGabor)).^2 + (yy-pos(2,kGabor)).^2)/2/sigma^2);
                s=cos(N.stim.sf*2*pi*ramp);
                cg=.25*g.*s;
                cg(abs(cg)<1/128)=0; % correct for bit-depth of monitor
                Cframe=Cframe+cg;
            end
            
            imagesc(xax, yax, Cframe, [-1 1]); colormap gray
            axis xy
            xlim([mnx mxx])
            ylim([mny mxy])
            hold on
            
            plot(0,0,'r+', 'MarkerFaceColor', 'r')
            t1XY=mode(N.stim.targ1XY);
            t2XY=mode(N.stim.targ2XY);
            plot(t1XY(1), t1XY(2), 'ro', 'MarkerFaceColor', 'r');
            plot(t2XY(1), t2XY(2), 'bo', 'MarkerFaceColor', 'b');
            
            grid on
            set(gca, 'gridcolor', 'w')
            
        end
       
        
        %--- plot neuron RF
        function h=plotMap(N, varargin)
            p=inputParser();
            p.addOptional('plotTask', false)
            p.addOptional('thetaRange', 360)
            p.parse(varargin{:});
            h=gca;
            if N.isMT
                if ~isempty(N.nreader.mtrfmap)
                    
                    map=N.nreader.mtrfmap;
                    
                    
%                     thetas=map.thetas;
                    
                    th=N.stim.theta;
                    thetas=map.thetas;
%                     if p.Results.plotTask
                        
%                     thetas=map.thetas-th;
                    if p.Results.thetaRange==360
                        
                    else
                        thetas(thetas>180)=thetas(thetas>180)-360;
                        thetas(thetas<-180)=thetas(thetas<-180)+360;
                    end
                    
                    rateMu=map.rateMu;
                    rateSd=map.rateStd;
                    [thetas, ix] = sort(thetas);
                    rateMu=rateMu(ix);
                    rateSd=rateSd(ix);
                    
                    vonmises=@(param,x)min(rateMu) +(max(rateMu)-min(rateMu))*exp(-param(2)*(1-cosd(x-param(1))));
                    param0=[circ_mean(thetas/180*pi, rateMu)/pi*180  circ_var(thetas/180*pi, rateMu)/pi*180  min(rateMu)];
                    evalc('param=lsqcurvefit(vonmises, param0, thetas, rateMu);');
                    
                    cmap=getMTLIPcolors;
                    if p.Results.thetaRange==360
                    plot(0:360, vonmises(param, 0:360), 'Color', cmap(1,:), 'Linewidth', 1);
                    xlim([0 360])
                    set(gca, 'Xtick', 0:90:360)
                    else
                    plot(-180:180, vonmises(param, -180:180), 'Color', cmap(1,:), 'Linewidth', 1);
                    xlim([-180 180])
                    set(gca, 'Xtick', -180:90:180)
                    end
                    hold on
                    
                    h=errorbar([thetas(end); thetas], [rateMu(end); rateMu],[rateSd(end); rateSd], '.k', 'MarkerSize', 5); hold on
                    plot([0 0], ylim, 'Color', .9*[1 1 1]);
                    plot(90*[1 1], ylim, 'Color', .9*[1 1 1]);
                    plot(-90*[1 1], ylim, 'Color', .9*[1 1 1]);
                    
%                     plot(-180:180, vonmises(param, -180:180), 'Color', get(h,'Color'))
                    
                    
                    xlabel('\Delta \theta')
                    ylabel('sp s^{-1}')
                    
                    
                elseif ~isempty(N.nreader.hyperflow)
                    map=N.nreader.hyperflow;
                    
                    
                    dy = repmat(sind(N.stim.theta), N.stim.nGabors,1);
                    dx = repmat(cosd(N.stim.theta), N.stim.nGabors, 1);
                    g = quiver(N.stim.gaborXY(:,1), N.stim.gaborXY(:,2), dx, dy, 'Color', 'r');
                    
                    set(g, 'AutoScaleFactor', .5)
                    %                     g = plot(N.stim.gaborXY(:,1), N.stim.gaborXY(:,2), 'o', 'MarkerFaceColor', 'r', 'Linestyle', 'none', 'MarkerSize', 10);
                    axis tight
                    hold on
                    nq = quiver(map.gridx(:), ...
                        map.gridy(:), ...
                        map.dx(:), ...
                        map.dy(:), 'k');
                    set(nq, 'AutoScaleFactor', .5);
                    axis off
                    
                    
                end
            else
                if ~isempty(N.nreader.delayedsaccades)
                    map=N.nreader.delayedsaccades;
                    cmap=getMTLIPcolors;
                    contourf(map.xax, map.yax, map.RF, [0 .5:.1:1], 'Color', cmap(2,:));
                    
%                     imagesc(map.xax, map.yax, map.RF);
                    colormap gray
                    hold on
                    xlim([-15 15]);
                    ylim([-12 12]);
                    set(gca, 'Xtick', [-20 -15 -10 -5 0 5 10 15 20], 'Ytick', [-10 -5 0 5 10])
                    set(gca, 'Color', [0 0 0])
                    %             rcg.plotTaskGeometry(stim, gca)
                    if p.Results.plotTask
                        plot(mode(N.stim.targ1XY(:,1)),mode(N.stim.targ1XY(:,2)), 'or')
                        plot(mode(N.stim.targ2XY(:,1)),mode(N.stim.targ2XY(:,2)), 'ob')
                        plot(mean(N.stim.gaborXY(:,1)), mean(N.stim.gaborXY(:,2)), 'og')
%                         axis off
                    end
                    plot([0 0], ylim, 'Color', .9*[1 1 1])
                    plot(xlim,[0 0], 'Color', .9*[1 1 1])
                end
                xlabel('space (deg)')
                ylabel('space (deg)')
                
            end
        end
        
        %--- does the neuron have mapping data?
        function b = hasMap(N)
            b = false;
            if N.isMT
                hmap = N.nreader.hyperflow;
                dmap = N.nreader.mtrfmap;
                if ~isempty(hmap) || ~isempty(dmap)
                    b = true;
                    return
                end
                
            elseif N.isLIP
                smap = N.nreader.delayedsaccades;
                b= ~isempty(smap);
                return
            end
        end 
        
        %--- does the neuron have simultaneously recorded MT neurons
        function b=hasMT(N)
            otherNeurons=setdiff(findFile(N.dirOld, N.exname), N.filename);
            brarea=cell(numel(otherNeurons), 1);
            for k=1:numel(otherNeurons)
                load(fullfile(N.dirOld, otherNeurons{k}), 'brainArea')
                brarea{k}=brainArea;
            end
            b=sum(strcmp(brarea, 'MT'));
        end
        
        %--- does the neuron have simultaneously recorded LIP neurons
        function b=hasLIP(N)
            otherNeurons=setdiff(findFile(N.dirOld, N.exname), N.filename);
            brarea=cell(numel(otherNeurons), 1);
            for k=1:numel(otherNeurons)
                load(fullfile(N.dirOld, otherNeurons{k}), 'brainArea')
                brarea{k}=brainArea;
            end
            b=sum(strcmp(brarea, 'LIP'));
        end
        
        function dpcell=getDprime(N)
            % optional arguments
%             'win', [-.2 2])
%             'binSize', .01)
%             'sort', 'net')
%             'smoothing', 0)
%             'aligningField', 'motionon')
            
            win = [-.2 2];
            bs  = .01;
            aligningField='motionon';
            ds=load(N.stim.filename, 'pulses', 'timing', 'targchosen', 'targcorrect', 'goodtrial');
            
            trialIndex=N.nreader.trialIndex;
            trialIndex=trialIndex(ds.goodtrial(trialIndex));
            ev=[ds.timing(trialIndex).(aligningField)]+[ds.timing(trialIndex).plxstart];
            trialIndex=trialIndex(~isnan(ev));
            ev=[ds.timing(trialIndex).(aligningField)]+[ds.timing(trialIndex).plxstart];
            
            fkern=[];
           
            [~,~,psthTime,~,r]=eventPsth(N.nreader.spikeTimes,ev,win,bs,fkern);
            
            assert(~any(any(isnan(r),2)), 'NaNs found')
                
            
            % stimulus for conditions
            pulses=mean(ds.pulses(trialIndex,:,:),3);          
            cho=ds.targchosen(trialIndex)==rcg.getTargPlus(ds);
            
            % get dprime during the first second of motion
            dpcell=arrayfun(@(x) dprime(mean(r(:,psthTime>0 & psthTime < x),2), sum(pulses,2)>0), 1.1);
%             dpcell=dprime(mean(r(:,psthTime>0 & psthTime < 1),2),sum(pulses,2)>0);
        end
        
        %% coherence sorted PSTH
        function S=coherencePSTH(N, varargin)
            % optional arguments
%             'win', [-.2 2])
%             'binSize', .01)
%             'sort', 'net')
%             'smoothing', 0)
%             'aligningField', 'motionon')
            
            p=inputParser();
            p.addOptional('win', [-.2 2])
            p.addOptional('binSize', .01)
            p.addOptional('sort', 'net')
            p.addOptional('smoothing', 0)
            p.addOptional('aligningField', 'motionon')
            p.addOptional('FixTarg', false)
            p.addOptional('targPlus', 1)
            p.addOptional('trialIndex', N.nreader.trialIndex);
            p.parse(varargin{:})
            
            win = p.Results.win;
            bs  = p.Results.binSize;
            aligningField=p.Results.aligningField;
            ds=load(N.stim.filename, 'pulses', 'timing', 'targchosen', 'targcorrect', 'goodtrial');
            
            trialIndex=p.Results.trialIndex;
            trialIndex=trialIndex(ds.goodtrial(trialIndex));
            ev=[ds.timing(trialIndex).(aligningField)]+[ds.timing(trialIndex).plxstart];
            trialIndex=trialIndex(~isnan(ev));
            ev=[ds.timing(trialIndex).(aligningField)]+[ds.timing(trialIndex).plxstart];
            xx=0:bs:10*bs;
            sig= p.Results.smoothing;
            if sig >0
%                 fkern=exp(- (xx.^2)/(2*sig^2)); 
                fkern=ones(sig,1);
                fkern=fkern/sum(fkern);
                    
            else
                fkern=[];
            end
            
            [~,~,psthTime,~,r]=eventPsth(N.nreader.spikeTimes,ev,win,bs,fkern);
            nPsthBins=numel(psthTime);
            assert(~any(any(isnan(r),2)), 'NaNs found')
                
            
            % stimulus for conditions
            pulses=mean(ds.pulses(trialIndex,:,:),3);
            pulses=pulses./std(pulses(:));
            
            if p.Results.targPlus==1
                cho=ds.targchosen(trialIndex)==rcg.getTargPlus(ds);
            else
                cho=ds.targchosen(trialIndex)==1;
            end
            correct=ds.targchosen(trialIndex)==ds.targcorrect(trialIndex);
            % get dprime during the first second of motion
            if strcmp(aligningField, 'motionon')
                dpcell=dprime(mean(r(:,psthTime>0 & psthTime < 1),2),cho);
            else
                dpcell=dprime(mean(r(:,psthTime>-.4 & psthTime < 0),2),cho);
            end
            
            if dpcell<0 && ~p.Results.FixTarg
                pulses=-pulses;
                cho=~cho;
%                 dpcell=dprime(mean(r(:,psthTime>0 & psthTime < 1),2),cho);
            end
            
            
            if ischar(p.Results.sort)
                coh=mean(pulses,2);
                abscoh=abs(coh);
                abscoh(abscoh==0)=.01;
                scoh=sign(cho-.5).*abscoh; % aligned to choice now
%                 binEdges = [-5 -1 -.5 -.25 0 .25 .5 1 5];
                binEdges = [-5 -1 -.25 0 .25 1 5];
                binCenters=binEdges(1:end-1)+diff(binEdges)/2;
                nCohs=numel(binEdges)-1;
                binId=cell2mat(arrayfun(@(x,y) scoh>x & scoh<y, binEdges(1:end-1), binEdges(2:end), 'UniformOutput', false));
            else
                kPulse=p.Results.sort;
                coh=pulses(:,kPulse);
                binCenters=unique(coh);
                nCohs=numel(binCenters);
                binId=cell2mat(arrayfun(@(x) coh==x, binCenters', 'UniformOutput', false));
            end
                
            
            
            
            
            psthCoh=zeros(nPsthBins, nCohs);
            psthCohCorr=zeros(nPsthBins,nCohs);
            
            for k=1:nCohs
                psthCoh(:,k)=mean(r(binId(:,k),:))/bs;
                psthCohCorr(:,k)=mean(r(binId(:,k)&correct,:))/bs;
            end
            
            S.psthTime=psthTime;
            S.psthCoh=psthCoh;
            S.psthCohCorr=psthCohCorr;
            S.psthCho = [mean(r(cho,:))' mean(r(~cho,:))']/bs;
            S.Cohs=binCenters;
            S.cohTrials=sum(binId);
            
        end

function g=getModelFit(N,modelDir, modelName, modelTag)
%g=getModelFit(N,modelDir, modelName, modelTag)
if nargin==1
    help neuro/getModelFit
    if nargout>0
        g=[];
    end
    return
end
neuronFits=findFile(modelDir, N.getName);
modelFits=findFile(neuronFits, modelName);
modelComparison=findFile(modelFits, 'modelComparison');
if ~isempty(modelComparison)
    modelFits(strcmp(modelFits, modelComparison))=[];
end
end

        

        
  
        
        %% WHITENED PULSE TRIGGERED AVERAGE
        function S=whitenedPTA(N)
            
            pulses=sum(N.stim.pulses,3); % single trial pulse sequences
%             pulses=zscore(pulses); % normalize to compare across datasets
            ev=N.stim.timing(:).motionon'+N.stim.timing(:).plxstart';
            cho=N.stim.targchosen==N.nreader.targPref;
            good=intersect(find(N.stim.goodtrial & ~isnan(ev)), N.nreader.trialIndex);
            if N.stim.targPlus~=N.nreader.targPref
                pulses=-pulses;
            end
            
            if strcmp(N.nreader.brainArea, 'LIP')
                binSize=20e-3;
%                 nkt=100;
                nkt=50;
            else
                binSize=10e-3;
                nkt=55;
            end
%             win=[-.1 3.5];
            win=[-.5 2.5];
            sig=.05;
%             fkernel=exp(-(-10:10).^2/(2*sig/binSize).^2); fkernel=fkernel/sum(fkernel);
%             fkernel=exp(-(1:100)/(sig/binSize)); fkernel=fkernel/sum(fkernel);
            fkernel=[];
            [~,~,bc,~,rtrial]=eventPsth(N.nreader.spikeTimes, ev(good), win, binSize, fkernel);
            rfilt=rtrial;
            rtc = rfilt;
            rtc(cho(good),:)  = repmat(mean(rfilt(cho(good),:)), sum(cho(good)),1);
            rtc(~cho(good),:) = repmat(mean(rfilt(~cho(good),:)), sum(~cho(good)),1);
            
            [~,i0]=min(bc.^2);
%             i0=0;
            ptimes=N.stim.timing(good(1)).pulses-N.stim.timing(good(1)).motionon;
            ptimes=ceil(ptimes/binSize)+i0;
            ptaTime=bsxfun(@plus, ((1:nkt)-i0)'*binSize, ptimes'*binSize);
            rPTA=pulseSTA(pulses(good,:),rfilt-rtc, ptimes, nkt, 1);
            rPTACho=pulseSTA(pulses(good,:),rtrial, ptimes, nkt, 1);
            
            ptaPulse=rPTA/binSize; % correct units (spikes/sec/coh)
            ptaCho=rPTACho/binSize;
            
            S.ptaTime=ptaTime;
            S.ptaPulse=ptaPulse;
            S.ptaCho=ptaCho;
            % smooth with half-gaussian kernel
            fkernel=exp(-(1:100).^2/(2*sig/binSize)^2); fkernel=fkernel/sum(fkernel);
            S.ptaPulseG=filter(fkernel,1,ptaPulse);
            S.ptaChoG=filter(fkernel,1,ptaCho);
            % smooth with exponential kernel
            fkernel=exp(-(1:100)/(sig/binSize)); fkernel=fkernel/sum(fkernel);
            S.ptaPulseE=filter(fkernel,1,ptaPulse);
            S.ptaChoE=filter(fkernel,1,ptaCho);
            S.snr=var(rPTACho(6:15,1))/var(rPTACho(1:5,1));
            S.ptaTrials=size(rtrial,1); % number of trials
        end
        
        
                %% WHITENED PULSE TRIGGERED AVERAGE Single Pulse Version
        function S=whitenedPTAsingle(N, type)
            
            if nargin<2
                type='none';
            end
            pulses=sum(N.stim.pulses,3); % single trial pulse sequences
%             pulses=zscore(pulses); % normalize to compare across datasets
            ev=N.stim.timing(:).motionon'+N.stim.timing(:).plxstart';
            cho=N.stim.targchosen==N.nreader.targPref;
            good=intersect(find(N.stim.goodtrial & ~isnan(ev)), N.nreader.trialIndex);
            if N.stim.targPlus~=N.nreader.targPref
                pulses=-pulses;
            end
            
            if strcmp(N.nreader.brainArea, 'LIP')
                binSize=20e-3;
%                 nkt=100;
                nkt=50;
            else
                binSize=10e-3;
                nkt=100;
            end
%             win=[-.1 3.5];
            win=[-.1 2.5];
            sig=.05;
%             fkernel=exp(-(-10:10).^2/(2*sig/binSize).^2); fkernel=fkernel/sum(fkernel);
%             fkernel=exp(-(1:100)/(sig/binSize)); fkernel=fkernel/sum(fkernel);
            fkernel=[];
            [~,~,bc,~,rtrial]=eventPsth(N.nreader.spikeTimes, ev(good), win, binSize, fkernel);
            rfilt=rtrial;
            
            
            [~,i0]=min(bc.^2);
%             i0=0;
            ptimes=N.stim.timing(good(1)).pulses-N.stim.timing(good(1)).motionon;
            ptimes=ceil(ptimes/binSize)+1;
            ptaTime=((1:nkt)-i0)'*binSize;
            
            switch type
                case 'none'
                    r=rfilt;
                case 'choice'
                    rtc = rfilt;
                    rtc(cho(good),:)  = repmat(mean(rfilt(cho(good),:)), sum(cho(good)),1);
                    rtc(~cho(good),:) = repmat(mean(rfilt(~cho(good),:)), sum(~cho(good)),1);
                    r=rfilt-rtc;
                case 'onset'
                    r0=ones(size(rfilt,1),1)*mean(rfilt);
                    r=rfilt-r0;
                otherwise
                    r=rfilt;
            end
                   
            r=r/binSize;
            rPTA=pulseSTASingle(pulses(good,:),r, ptimes, nkt);
            
            
            ptaPulse=rPTA; %/binSize; % correct units (spikes/sec/coh)
            
            S.ptaTime=ptaTime;
            S.ptaPulse=ptaPulse;
            % smooth with half-gaussian kernel
            fkernel=exp(-(1:100).^2/(2*sig/binSize)^2); fkernel=fkernel/sum(fkernel);
            S.ptaPulseG=filter(fkernel,1,ptaPulse);
            % smooth with exponential kernel
            fkernel=exp(-(1:100)/(sig/binSize)); fkernel=fkernel/sum(fkernel);
            S.ptaPulseE=filter(fkernel,1,ptaPulse);
            S.ptaTrials=size(rfilt,1); % number of trials
        end
        
        %% check the neuron overall
        function checkPlot(N)
            clf
            set(gcf, 'DefaultAxesColorOrder', hot(12))
            subplot(3,3,1)
            N.plotMap
            xlabel('')
            ylabel('')
            title('')
            subplot(3,3,2)
            pulses=mean(N.stim.pulses,3);
            ev=N.stim.timing(:).motionon'+N.stim.timing(:).plxstart';
            cho=N.stim.targchosen==N.nreader.targPref;
            good=intersect(find(N.stim.goodtrial & ~isnan(ev)), N.nreader.trialIndex);
            
            if N.stim.targPlus~=N.nreader.targPref
                pulses=-pulses;
            end
            
            if strcmp(N.nreader.brainArea, 'LIP')
                binSize=20e-3;
            else
                binSize=10e-3;
            end
            win=[-.5 2.5];
            sig=.025;
            fkernel=exp(-(-1:100).^2/(2*sig/binSize).^2); fkernel=fkernel/sum(fkernel);
            [~,~,bc,~,rtrial]=eventPsth(N.nreader.spikeTimes, ev(good), win, binSize);
            
            rfilt=filter(fkernel,1,rtrial')';
            
            m1=mean(rfilt(cho(good),:))/binSize;
            s1=std(rfilt(cho(good),:))/sqrt(sum(cho(good)))/binSize;
            m2=mean(rfilt(~cho(good),:))/binSize;
            s2=std(rfilt(~cho(good),:))/sqrt(sum(~cho(good)))/binSize;
            
%             rtc = rfilt;
%             rtc(cho(good),:)  = repmat(mean(rfilt(cho(good),:)), sum(cho(good)),1);
%             rtc(~cho(good),:) = repmat(mean(rfilt(~cho(good),:)), sum(~cho(good)),1);
            
            plot(bc, m1,'-', bc,m1+s1, '--', bc, m1-s1, '--', 'Color', [0    0.4470    0.7410]); hold on
            plot(bc, m2,'-', bc,m2+s2, '--', bc, m2-s2, '--', 'Color', [0.8500    0.3250    0.0980])
            xlim([-.2 1.5])
            
%             [~,i0]=min(bc.^2);
%             ptimes=N.stim.timing(good(1)).pulses-N.stim.timing(good(1)).motionon;
%             ptimes=ceil(ptimes/binSize)+i0;
            
            
%             nkt=50;
%             rPTA=pulseSTA(pulses(good,:),rfilt-rtc, ptimes, nkt, 1);
            P=N.whitenedPTA();
            subplot(3,3,3)
            plot(P.ptaTime, P.ptaPulseG)
            
            subplot(3,3,4:9)
            imagesc(bc,1:sum(cho(good)),-rtrial(cho(good),:)); colormap gray
            hold on
            imagesc(bc,sum(cho(good)):(sum(cho(good))+sum(~cho(good))),-rtrial(~cho(good),:)); colormap gray
            ylim([1 size(rtrial,1)])
            suptitle([N.nreader.exname '  ' num2str(N.nreader.id)])
            sciencestagram(gcf, 10, [8 11])
            figname = fullfile(getpref('mtlip', 'figures'), 'single_neuron_check', sprintf('%s_%02.0f.pdf', N.nreader.exname, N.nreader.id));
            tic
            print(gcf, '-dpdf', '-opengl', figname)
            toc
        end
        
        %% neuron's name for glm model fits
        function name=getName(N, includeExperiment)
            if nargin==1
                includeExperiment=true;
            end
            if includeExperiment
                name=sprintf('%s%sneuron%02.0fch%02.0f', N.nreader.exname, N.nreader.brainArea, N.nreader.id, N.nreader.channel);
            else
                name=sprintf('%sneuron%02.0fch%02.0f', N.nreader.brainArea, N.nreader.id, N.nreader.channel);
            end
        end
        
    end % public methods


end


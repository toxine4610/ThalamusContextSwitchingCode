%% parallelized for fast calculations
clear all;
addpath( genpath('ReliabTools\') );

load('C:\Users\halassalab\Dropbox\Rajeev\ForContextSwitchProject\AllMice_RuleSwitched_1002_SpkTimes.mat')

%%
delay = [0, 0.5];
baseline = [-1, 0];
ModulationIndex = @(x,y) (x-y)./(x+y);
MDColor = [255,102,102]./255;

%% set timing information for spike Rate calculations
pre = 0.5;
post = 1.5;

binNum = 120;
xlim_start = -0.5;
xlim_end = 1.5;


% plotting parameters
bin = 0.01;
filtWidth = 0.05;


%% sort into sessions

Type = 'PFC';
UnitNumbersPFC = [Spfc.unitNbr];
StartPointsPFC = find(UnitNumbersPFC == 1);

for i = 2:length(StartPointsPFC)+1
    if i <= length(StartPointsPFC)
        SessListPFC{i-1} = UnitNumbersPFC( StartPointsPFC(i-1) :  StartPointsPFC(i)-1 );
    else
        SessListPFC{i-1} = UnitNumbersPFC( StartPointsPFC(i-1): end );
    end;
    
end;

UnitNumbersMD = [Smd.unitNbr];
D = diff( UnitNumbersMD );
NZ = [0, find( D ~= 1)];
for i = 2:length(NZ)+1
    if i <= length(NZ)
        SessListMD{i-1} = UnitNumbersMD( NZ(i-1)+1 :  NZ(i) );
    else
        SessListMD{i-1}= UnitNumbersMD( NZ(i-1)+1 :  end );
    end;
end;

%% PFC

clear Consistency CCArray R1CMI R2CMI R1C1_RelTs R1C2_RelTs R2C1_RelTs R2C2_RelTs R1C1_maxRel R1C2_maxRel R2C1_maxRel R2C2_maxRel

% ppm = ParforProgMon('Processing PFC  ',length(UnitNumbersPFC));

%%
parfor s = 1:length(SessListPFC)
    units_this_session = SessListPFC{s} +  (StartPointsPFC(s)-1);
    ct=  0;
    
    for u = units_this_session
        
%         ppm.increment();
        fprintf('Processing Sess: %d Unit %d of %d.....', s, u, length(units_this_session));
        ct= ct+1;
        
        %%
        
        % ====  Extract spike times.....................................
        
        R1C1 = Spfc(u).SpikeTimes_R1C1; %1
        R1C2 = Spfc(u).SpikeTimes_R1C2; %2
        
        R2C1 = Spfc(u).SpikeTimes_R2C1; %3
        R2C2 = Spfc(u).SpikeTimes_R2C2; %4
        
        % ====  Compute Reliability TS ....................................
        [R1C1_weightedCC, R1C1_weightedCC_scaled] = computeRawReliability(R1C1);
        [R1C2_weightedCC, R1C2_weightedCC_scaled] = computeRawReliability(R1C2);
        
        [R2C1_weightedCC, R2C1_weightedCC_scaled] = computeRawReliability(R2C1);
        [R2C2_weightedCC, R2C2_weightedCC_scaled] = computeRawReliability(R2C2);
        
    
        pR1C1 = determineSignificantPeaks(R1C1);
        pR1C2 = determineSignificantPeaks(R1C2);
        pR2C1 = determineSignificantPeaks(R2C1);
        pR2C2 = determineSignificantPeaks(R2C2);
       
        P{s}(ct,:) = [pR1C1,pR1C2,pR2C1,pR2C2];
        
        %% Determine consistency
        % Rule consistent peaks
        % Rule - 1
        [R1_maxCC, R1_maxLag, R1_Consistent] = computeSimilarityPeaks( R1C1_weightedCC_scaled, R1C2_weightedCC_scaled );
        % Rule - 2
        [R2_maxCC, R2_maxLag, R2_Consistent] = computeSimilarityPeaks( R2C1_weightedCC_scaled, R2C2_weightedCC_scaled );
        
        % Context consistency
        % Context - 1
        [C1_maxCC, C1_maxLag, C1_Consistent] = computeSimilarityPeaks( R1C1_weightedCC_scaled, R2C1_weightedCC_scaled );
        % Context -2
        [C2_maxCC, C2_maxLag, C2_Consistent] = computeSimilarityPeaks( R1C2_weightedCC_scaled, R2C2_weightedCC_scaled );
        
        % Cue consistency
        % Cue - 1
        [Q1_maxCC, Q1_maxLag, Q1_Consistent] = computeSimilarityPeaks( R1C1_weightedCC_scaled, R2C2_weightedCC_scaled );
        % Cue -2
        [Q2_maxCC, Q2_maxLag, Q2_Consistent] = computeSimilarityPeaks( R2C1_weightedCC_scaled, R1C2_weightedCC_scaled );
        
        Consistency{s}(ct,:) = [R1_Consistent, R2_Consistent, C1_Consistent, C2_Consistent, Q1_Consistent, Q2_Consistent];
        CCArray{s}{ct}       = [R1_maxCC, R2_maxCC; C1_maxCC, C2_maxCC; Q1_maxCC, Q2_maxCC];
        
        %%
        [maxRelR1C1, maxRelR1C1_time] = max(R1C1_weightedCC);
        [maxRelR1C2, maxRelR2C1_time] = max(R1C2_weightedCC);
        
        [maxRelR2C1, maxRelR1C2_time] = max(R2C1_weightedCC);
        [maxRelR2C2, maxRelR2C2_time] = max(R2C2_weightedCC);
        
        R1CMI{s}(ct)  = ModulationIndex( maxRelR1C1, maxRelR1C2 );
        R2CMI{s}(ct)  = ModulationIndex( maxRelR2C1, maxRelR2C2 );
        
        
        R1C1_RelTs{s}{ct}  = R1C1_weightedCC;
        R1C2_RelTs{s}{ct}  = R1C2_weightedCC;
        
        R2C1_RelTs{s}{ct}  = R2C1_weightedCC;
        R2C2_RelTs{s}{ct}  = R2C2_weightedCC;
        
        R1C1_maxRel{s}(ct) =  maxRelR1C1;
        R2C1_maxRel{s}(ct) =  maxRelR2C1;
        R1C2_maxRel{s}(ct) =  maxRelR1C2;
        R2C2_maxRel{s}(ct) =  maxRelR2C2;
        %         clear R1C1_weightedCC R1C1_weightedCC_scaled R1C2_weightedCC R1C2_weightedCC_scaled ...
        %               R2C1_weightedCC R2C1_weightedCC_scaled R2C2_weightedCC R2C2_weightedCC_scaled
        
        fprintf('...Done\n');
    end;
end;

PFC.P            = P;
PFC.Consistency  = Consistency;
PFC.CCArray      = CCArray;
PFC.R1CMI        =   R1CMI;
PFC.R2CMI        =   R2CMI;

PFC.R1C1_RelTs   =  R1C1_RelTs;
PFC.R1C2_RelTs   =  R1C2_RelTs;
PFC.R1C1_maxRel  =  R1C1_maxRel;
PFC.R2C1_maxRel  =  R2C1_maxRel;

PFC.R2C1_RelTs   =  R2C1_RelTs;
PFC.R2C2_RelTs   =  R2C2_RelTs;
PFC.R1C2_maxRel  =  R1C2_maxRel;
PFC.R2C2_maxRel  =  R2C2_maxRel;

clear Consistency CCArray R1CMI R2CMI R1C1_RelTs R1C2_RelTs R2C1_RelTs R2C2_RelTs R1C1_maxRel R1C2_maxRel R2C1_maxRel R2C2_maxRel


%% MD

clear Consistency  P CCArray R1CMI R2CMI R1C1_RelTs R1C2_RelTs R2C1_RelTs R2C2_RelTs R1C1_maxRel R1C2_maxRel R2C1_maxRel R2C2_maxRel


% ppm = ParforProgMon('Processing MD  ',length(UnitNumbersMD));


parfor s = 1:length(SessListMD)
    
    
    units_this_session = SessListMD{s}-( (SessListMD{s}(1)-1) ) + NZ(s);
    ct=  0;
    for u = units_this_session
        %         clear R1C1 R1C2 R2C1 R2C2
%         ppm.increment();
        fprintf('Processing MD:  Sess: %d Unit %d of %d.....', s, u, length(units_this_session));
        ct= ct+1;
        
        %%
        
        % ====  Extract spike times.....................................
        
        R1C1 = Smd(u).SpikeTimes_R1C1; %1
        R1C2 = Smd(u).SpikeTimes_R1C2; %2
        
        R2C1 = Smd(u).SpikeTimes_R2C1; %3
        R2C2 = Smd(u).SpikeTimes_R2C2; %4
        
        % ====  Compute Reliability TS ....................................
        [R1C1_weightedCC, R1C1_weightedCC_scaled] = computeRawReliability(R1C1);
        [R1C2_weightedCC, R1C2_weightedCC_scaled] = computeRawReliability(R1C2);
        
        [R2C1_weightedCC, R2C1_weightedCC_scaled] = computeRawReliability(R2C1);
        [R2C2_weightedCC, R2C2_weightedCC_scaled] = computeRawReliability(R2C2);
        
         pR1C1 = determineSignificantPeaks(R1C1);
        pR1C2 = determineSignificantPeaks(R1C2);
        pR2C1 = determineSignificantPeaks(R2C1);
        pR2C2 = determineSignificantPeaks(R2C2);
       
        P{s}(ct,:) = [pR1C1,pR1C2,pR2C1,pR2C2];
        
        %% Determine consistency
        % Rule consistent peaks
        % Rule - 1
        [R1_maxCC, R1_maxLag, R1_Consistent] = computeSimilarityPeaks( R1C1_weightedCC_scaled, R1C2_weightedCC_scaled );
        % Rule - 2
        [R2_maxCC, R2_maxLag, R2_Consistent] = computeSimilarityPeaks( R2C1_weightedCC_scaled, R2C2_weightedCC_scaled );
        
        % Context consistency
        % Context - 1
        [C1_maxCC, C1_maxLag, C1_Consistent] = computeSimilarityPeaks( R1C1_weightedCC_scaled, R2C1_weightedCC_scaled );
        % Context -2
        [C2_maxCC, C2_maxLag, C2_Consistent] = computeSimilarityPeaks( R1C2_weightedCC_scaled, R2C2_weightedCC_scaled );
        
        % Cue consistency
        % Cue - 1 (blue)
        [Q1_maxCC, Q1_maxLag, Q1_Consistent] = computeSimilarityPeaks( R1C1_weightedCC_scaled, R2C2_weightedCC_scaled );
        % Cue -2 (brown)
        [Q2_maxCC, Q2_maxLag, Q2_Consistent] = computeSimilarityPeaks( R2C1_weightedCC_scaled, R1C2_weightedCC_scaled );
        
        Consistency{s}(ct,:) = [R1_Consistent, R2_Consistent, C1_Consistent, C2_Consistent, Q1_Consistent, Q2_Consistent];
        CCArray{s}{ct}       = [R1_maxCC, R2_maxCC; C1_maxCC, C2_maxCC; Q1_maxCC, Q2_maxCC];
        
        %%
        [maxRelR1C1, maxRelR1C1_time] = max(R1C1_weightedCC);
        [maxRelR1C2, maxRelR2C1_time] = max(R1C2_weightedCC);
        
        [maxRelR2C1, maxRelR1C2_time] = max(R2C1_weightedCC);
        [maxRelR2C2, maxRelR2C2_time] = max(R2C2_weightedCC);
        
        R1CMI{s}(ct)  = ModulationIndex( maxRelR1C1, maxRelR1C2 );
        R2CMI{s}(ct)  = ModulationIndex( maxRelR2C1, maxRelR2C2 );
        
        
        R1C1_RelTs{s}{ct}  = R1C1_weightedCC;
        R1C2_RelTs{s}{ct}  = R1C2_weightedCC;
        
        R2C1_RelTs{s}{ct}  = R2C1_weightedCC;
        R2C2_RelTs{s}{ct}  = R2C2_weightedCC;
        
        R1C1_maxRel{s}(ct) =  maxRelR1C1;
        R2C1_maxRel{s}(ct) =  maxRelR2C1;
        R1C2_maxRel{s}(ct) =  maxRelR1C2;
        R2C2_maxRel{s}(ct) =  maxRelR2C2;
        %         clear R1C1_weightedCC R1C1_weightedCC_scaled R1C2_weightedCC R1C2_weightedCC_scaled ...
        %               R2C1_weightedCC R2C1_weightedCC_scaled R2C2_weightedCC R2C2_weightedCC_scaled
        
        fprintf('...Done\n');
    end;
end;

MD.P            = P;
MD.Consistency  = Consistency;
MD.CCArray      = CCArray;
MD.R1CMI        =   R1CMI;
MD.R2CMI        =   R2CMI;

MD.R1C1_RelTs   =  R1C1_RelTs;
MD.R1C2_RelTs   =  R1C2_RelTs;
MD.R1C1_maxRel  =  R1C1_maxRel;
MD.R2C1_maxRel  =  R2C1_maxRel;

MD.R2C1_RelTs   =  R2C1_RelTs;
MD.R2C2_RelTs   =  R2C2_RelTs;
MD.R1C2_maxRel  =  R1C2_maxRel;
MD.R2C2_maxRel  =  R2C2_maxRel;

clear Consistency CCArray R1CMI R2CMI R1C1_RelTs R1C2_RelTs R2C1_RelTs R2C2_RelTs R1C1_maxRel R1C2_maxRel R2C1_maxRel R2C2_maxRel


% %%
% 
% 
% 
% PeakSimIndex = cell2mat( PFC.Consistency(1:3)' );
% 
% R1Const      = find( PeakSimIndex(:,1) == 1);
% R2Const      = find( PeakSimIndex(:,2) == 1);
% CommonMode   = intersect( R1Const, R2Const);
% R1Const      = setdiff( R1Const, CommonMode);
% R2Const      = setdiff( R2Const, CommonMode) ;
% R1R2Const    = CommonMode;
% 
% R1CMIpool   =  cell2mat(PFC.R1CMI );
% R2CMIpool   = cell2mat(PFC.R2CMI );
% 
% CMI_Rel_Rule1   = R1CMIpool( R1Const );
% CMI_Rel_Rule2   = R2CMIpool( R2Const );
% 
% 
% 
% 
% 
% %%
% Rate = load('Rate_CS.mat');
% CMI_Rate_Rule1 = cell2mat( Rate.PFC.MI_mu_12(1:3) );
% CMI_Rate_Rule2 = cell2mat( Rate.PFC.MI_mu_34(1:3) );
% 
% CMI_Rate_Rule1 = CMI_Rate_Rule1( R1Const );
% CMI_Rate_Rule2 = CMI_Rate_Rule2( R2Const );
% 
% 
% 
% 
% %%
% figure(1); set(gcf,'color','w');
% subplot(1,2,1);
% plot( CMI_Rate_Rule1, CMI_Rel_Rule1,'o','markersize',12, 'markerfacecolor','w','markeredgecolor','k'); hold on;
% axis square; box off;
% set(gca,'xtick',-1:0.2:1, 'ytick',-1:0.2:1)
% set(gca,'tickdir','out','fontsize',15); box off;
% line( [ 0,0], [-1,1],'color','k');
% line( [-1, 1], [0, 0 ],'color','k');
% plot( linspace(-1, 1, 20000), linspace(-1, 1, 20000),'color','k','linewidth',3)
% 
% 
% subplot(1,2,2);
% plot( CMI_Rate_Rule2, CMI_Rel_Rule2,'o','markersize',12, 'markerfacecolor','w','markeredgecolor','k'); hold on
% axis square; box off;
% set(gca,'tickdir','out','fontsize',15); box off;
% set(gca,'xtick',-1:0.2:1, 'ytick',-1:0.2:1)
% line( [ 0,0], [-1,1],'color','k');
% line( [-1, 1], [0, 0 ],'color','k');
% plot( linspace(-1, 1, 20000), linspace(-1, 1, 20000),'color','k','linewidth',3)


context = Z_C1(:,9);
dC      = diff(context);
SwitchTrial = find(dC~=0); % switch occurs after this trial

if length(SwitchTrial) == 1
    RangeC1 = 1:SwitchTrial;
    RangeC2 = SwitchTrial+1:size(Z_C1,1);
elseif length(SwitchTrial) == 2
    RangeC1 = 1:SwitchTrial(1);
    RangeC2 = SwitchTrial(1)+1 : SwitchTrial(2);
    RangeC1rep = SwitchTrial(2)+1: size(Z_C1,1);
end;


Z_context1 = Z_C1(RangeC1,:);
Z_context2 = Z_C1(RangeC2,:);

%%
ct1 = []; ct2 = [];ct3 = []; ct4 = [];
for t = 1:size(Z_C1,1)
    if Z_C1(t,1) == 1 && Z_C1(t,2) == 1
        ct1 = [ct1, t];
    elseif Z_C1(t,1) == 0 && Z_C1(t,2) == 1
        ct2 = [ct2, t];
        
    elseif Z_C1(t,1) == 1 && Z_C1(t,2) == 2
        ct3 = [ct3, t];
    elseif Z_C1(t,1) == 0 && Z_C1(t,2) == 2
        ct4 = [ct4, t];
    end;
end;

%% compile all spikes.....................................................

for i= 1:numel(Spfc)
    
    SpikeTimes = cell(size(Z_C1,1),1);
    
    trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 1;
    
    SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R1C1;
    
    trialSubset = Z_C1(:,1)==0 & Z_C1(:,2) == 1;
    
    SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R1C1_inc;
    
    trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 2;
    
    SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R2C1;
    
    trialSubset = Z_C1(:,1)==0 & Z_C1(:,2) == 2;
    
    SpikeTimes(trialSubset) = Spfc(i).SpikeTimes_R2C1_inc;
    
    PFC(i).SpikeTimes = SpikeTimes;
    [spikeRate, spikeMat, timeOut, indBase] = makeSpikeRates(PFC(i).SpikeTimes, [-0.1, 0.7] , 0.001, 0.03);
    PFC(i).muRate = spikeRate;
    PFC(i).Rate   = spikeMat;
    PFC(i).TRRate = nanmean( spikeMat(:, find(timeOut >= 0 & timeOut <= 0.6 )),2);
    PFC(i).Ratenorm = ( spikeRate - nanmin( spikeRate ) ) ./ ( nanmax( spikeRate ) - nanmin( spikeRate ) );
    PFC(i).TRRatenorm = ( PFC(i).TRRate - nanmin(PFC(i).TRRate ) ) ./ ( nanmax( PFC(i).TRRate ) - nanmin( PFC(i).TRRate ) );
    
    [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(PFC(i).SpikeTimes );
    R(i) = findFractionBlankTrials(SpikeCount, bins_for_plotting);
    
end;

goodPFC = R < 0.55;
PFC = PFC(goodPFC);



clear R;
for i= 1:numel(Smd)
    
    SpikeTimes = cell(size(Z_C1,1),1);
    
    trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 1 ;
    
    SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R1C1;
    
    trialSubset = Z_C1(:,1)==0 & Z_C1(:,2) == 1;
    
    SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R1C1_inc;
    
    trialSubset = Z_C1(:,1)==1 & Z_C1(:,2) == 2;
    
    SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R2C1;
    
    trialSubset = Z_C1(:,1)==0 & Z_C1(:,2) == 2;
    
    SpikeTimes(trialSubset) = Smd(i).SpikeTimes_R2C1_inc;
    
    MD(i).SpikeTimes = SpikeTimes;
    [spikeRate, spikeMat, timeOut, indBase] = makeSpikeRates(MD(i).SpikeTimes, [-0.1, 0.7] , 0.001, 0.03);
    MD(i).muRate = spikeRate;
    MD(i).Rate   = spikeMat;
    MD(i).TRRate = nanmean( spikeMat(:, find(timeOut >= 0 & timeOut <= 0.6 )),2);
    MD(i).Ratenorm = ( spikeRate - nanmin( spikeRate ) ) ./ ( nanmax( spikeRate ) - nanmin( spikeRate ) );
    MD(i).TRRatenorm = ( MD(i).TRRate - nanmin(MD(i).TRRate ) ) ./ ( nanmax( MD(i).TRRate ) - nanmin(MD(i).TRRate ) );
    
     [FF_SC, mean_SC, var_SC, bins_for_plotting, SpikeCount, ISI] = computeFanoFactor(MD(i).SpikeTimes );
    R(i) = findFractionBlankTrials(SpikeCount, bins_for_plotting);
end;

goodMD = R < 0.55;
MD = MD(goodMD);
%%

V = [PFC.TRRate];

win = 4; len = 10;
ct = 0;

% ============= correct laser off;
clear mu numC
for i = 1:win
    ct = ct+1;
    range = (SwitchTrial(1) + (i-(win-1))*len) :SwitchTrial(1) + (i-(win-2))*len
    ind = Z_C1(range, 1) == 1 & Z_C1(range, 10) == 0;
    
    foo    =   V(range(ind),:);
    if size(foo,1) > 1
        mu{ct} =  mean(V(range(ind),:),1)';
    elseif size(foo,1) == 1;
        mu{ct} =  V(range(ind),:)';
    end;
    
    numC(ct) = length( find(Z_C1(range,1)==1) )./length(range);
    
end;
   
figure(1);
subplot(1,4,1); 
title('Correct, Laser Off')
y  = cellfun(@(x) nanmean(x), mu );
yyaxis left;
plot( 1:ct, y,'color','r'); hold on;
yyaxis right;
plot( 1:ct, numC,'-o','color','b','linewidth',3); hold on;
axis square; box off; set(gca,'tickdir','out','fontsize',16);
xlim([0.5, 4.5]);



%%
% ============= correct laser on;
ct = 0;
clear mu numC
for i = 1:win
    ct = ct+1;
    range = (SwitchTrial(2) + (i-(win-1))*len) :SwitchTrial(2) + (i-(win-2))*len
    ind = Z_C1(range, 1) == 1 & Z_C1(range, 10) == 1;
    
    foo    =   V(range(ind),:);
    if size(foo,1) > 1
        mu{ct} =  mean(V(range(ind),:),1)';
    elseif size(foo,1) == 1;
        mu{ct} =  V(range(ind),:)';
    end;
    
    numC(ct) = length( find(Z_C1(range,1)==1))./length(range);
    
end;
   
figure(1);
subplot(1,4,2);
title('Correct, Laser On')
y  = cellfun(@(x) nanmean(x), mu );
yyaxis left;
plot( 1:ct, y,'color','r'); hold on;
yyaxis right;
plot( 1:ct, numC,'-o','color','b','linewidth',3); hold on;
axis square; box off; set(gca,'tickdir','out','fontsize',16);
xlim([0.5, 4.5]);


% ============= wrong laser off;

ct = 0;
clear mu numC
for i = 1:win
    ct = ct+1;
    range = (SwitchTrial + (i-(win-1))*len) :SwitchTrial + (i-(win-2))*len
    ind = Z_C1(range, 1) == 0 & Z_C1(range, 10) == 0;
    
    foo    =   V(range(ind),:);
    if size(foo,1) > 1
        mu{ct} =  mean(V(range(ind),:),1)';
    elseif size(foo,1) == 1;
        mu{ct} =  V(range(ind),:)';
    end;
    
    numC(ct) = length( find(Z_C1(range,1)==1))./length(range);
    
end;
   
figure(1);
subplot(1,4,3);
title('Wrong, Laser Off')
y  = cellfun(@(x) nanmean(x), mu );
yyaxis left;
plot( 1:ct, y,'color','r'); hold on;
yyaxis right;
plot( 1:ct, numC,'-o','color','b','linewidth',3); hold on;
axis square; box off; set(gca,'tickdir','out','fontsize',16);
xlim([0.5, 4.5]);

% ============= wrong laser on;

ct = 0;
for i = 1:win
    ct = ct+1;
    range = (SwitchTrial + (i-(win-1))*len) :SwitchTrial + (i-(win-2))*len
    ind = Z_C1(range, 1) == 0 & Z_C1(range, 10) == 1;
    
    foo    =   V(range(ind),:);
    if size(foo,1) > 1
        mu{ct} =  mean(V(range(ind),:),1)';
    elseif size(foo,1) == 1;
        mu{ct} =  V(range(ind),:)';
    end;
    
    numC(ct) = length( find(Z_C1(range,1)==1))./length(range);
    
end;
   
figure(1);
subplot(1,4,4);
title('Wrong, Laser On')
y  = cellfun(@(x) nanmean(x), mu );
yyaxis left;
plot( 1:win, y,'color','r'); hold on;
yyaxis right;
plot( 1:win, numC,'-o','color','b','linewidth',3); hold on;
axis square; box off; set(gca,'tickdir','out','fontsize',16);
xlim([0.5, 4.5]);



%%

fisher  = @(r) 0.5.*log( (1+r)./(1-r) );


win = 4; len = 10;
P = nchoosek(1:numel(MD),2);
clear CC
for i = 1:win
    
    range = (SwitchTrial(1) + (i-(win-1))*10) :SwitchTrial(1) + (i-(win-2))*10
    ind = Z_C1(range, 1) == 0 ;
    for j = 1:size(P,1)
        clear cell1 cell2
        cell1 = MD(P(j,1)).Rate(range(ind),find(timeOut >= 0 & timeOut <= 0.6 ));
        cell2 = MD(P(j,2)).Rate(range(ind),find(timeOut >= 0 & timeOut <= 0.6 ));
    
        CC(i,j) = corr( nanmean(cell1)', nanmean(cell2)');    
    end;
end;

figure(2);
boxplot(CC'); hold on;
plotSpread(CC');
axis square; box off; set(gca,'tickdir','out','fontsize',16);

figure(3);
yyaxis left;
plot( 1:win, nanmedian(CC,2),'color','g'); hold on;
yyaxis right;
plot( 1:win, numC,'-o','color','b','linewidth',3); hold on;
axis square; box off; set(gca,'tickdir','out','fontsize',16);
xlim([0.5, 4.5]);

%%

for i = 1:numel(MD)
    this_cell = MD(i).Rate;
    ind = Z_C1(RangeC1, 1) == 1 ;
    C1Ave     = nanmean( this_cell( RangeC1(ind),:),1);
    
    ind = Z_C1(RangeC1rep, 1) == 1 ;
    C2Ave     = nanmean( this_cell( RangeC1rep(ind),:),1);
    for j = RangeC1
        Rel_C1(i,j) = corr( this_cell(j,:)', C1Ave' );
    end;
    for j = RangeC1rep
        Rel_C2(i,j-RangeC2(1)+1) = corr( this_cell(j,:)', C2Ave');
    end;
end


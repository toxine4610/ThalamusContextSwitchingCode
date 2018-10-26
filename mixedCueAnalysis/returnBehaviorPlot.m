
function returnBehaviorPlot( Z_C1, d)


context = Z_C1(:,9);
dC      = diff(context);
SwitchTrial = find(dC~=0); % switch occurs after this trial

if length(SwitchTrial) == 1
    RangeC1 = 1:SwitchTrial;
    RangeC2 = SwitchTrial+1:size(Z_C1,1);
    
    Z_context1 = Z_C1(RangeC1,:);
    Z_context2 = Z_C1(RangeC2,:);
    Z_context1rep = [];
elseif length(SwitchTrial) == 2
    RangeC1 = 1:SwitchTrial(1);
    RangeC2 = SwitchTrial(1)+1 : SwitchTrial(2);
    RangeC1rep = SwitchTrial(2)+1: size(Z_C1,1);
    Z_context1 = Z_C1(RangeC1,:);
    Z_context2 = Z_C1(RangeC2,:);
    Z_context1rep = Z_C1(RangeC1rep,:);
else
    Z_context1 = Z_C1;
end;

%%

w = gausswin(5);
w = w/sum(w);

perf1 =conv(Z_context1(:,1),w,'same');
perf2 =conv(Z_context2(:,1),w,'same');
perf3 =conv(Z_context1rep(:,1),w,'same');

P = [smooth(perf1) ; smooth(perf2) ; smooth(perf3) ];

figure(d); set(gcf,'color','w');
plot( 1:length(P), P, 'color','k','linewidth',3 ); hold on
for i = 1:length(SwitchTrial)
    line( [SwitchTrial(i), SwitchTrial(i)], [0,1],'color','r' );
end;

line([0,200],[0.5, 0.5],'color','r')
axis normal; box off; set(gca,'tickdir','out','fontsize',16);

%%

[ExpH_C1, rup_C1, rlo_C1] = computeExpectation(Z_context1);
[ExpH_C2, rup_C2, rlo_C2] = computeExpectation(Z_context2);
[ExpH_C3, rup_C3, rlo_C3] = computeExpectation(Z_context1rep);
[ExpH_All, rup, rlo] = computeExpectation(Z_C1);

figure(d);
plot( 1:length(ExpH_C1), (ExpH_C1), 'linewidth',3); hold on
plot( linspace(SwitchTrial(1), SwitchTrial(2), length(ExpH_C2)), (ExpH_C2),'linewidth',3 );
plot( linspace(SwitchTrial(2), max(size(Z_C1)), length(ExpH_C3)), (ExpH_C3),'linewidth',3 );
plot( 1:size(Z_C1,1), smooth(ExpH_All),'r','linewidth',3);

drawnow;

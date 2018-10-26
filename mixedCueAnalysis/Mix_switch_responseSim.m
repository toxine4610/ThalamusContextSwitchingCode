

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
Z_context1rep = Z_C1(RangeC1rep,:);


%%
[PFC, MD, goodPFC, goodMD] = cleanData(Spfc, Smd, Z_C1);

Spfc = Spfc(goodPFC);
Smd = Smd(goodMD);

%%

[Spfc_mix_C1, Smd_mix_C1, ZVis, ZAud] = packageData(Z_context1, Smd, Spfc);

[Spfc_mix_C1rep, Smd_mix_C1rep, ZVis, ZAud] = packageData(Z_context1rep, Smd, Spfc);

[Spfc_mix_C2, Smd_mix_C2, ZVis, ZAud] = packageData(Z_context2, Smd, Spfc);


C1    = extractVariables_PFC( Spfc_mix_C1 );
C1rep = extractVariables_PFC( Spfc_mix_C1rep );
C2    = extractVariables_PFC( Spfc_mix_C2 );


C1_MD    = extractVariables_MD( Smd_mix_C1 );
C1rep_MD = extractVariables_MD( Smd_mix_C1rep );
C2_MD    = extractVariables_MD( Smd_mix_C2 );


%%

CMI   = @(X,Y) (X-Y)./(X+Y);

cmi_block1 = CMI( C1.PFC_muFR_CFD_C1, C2.PFC_muFR_CFD_C1 );
cmi_block2 = CMI( C1rep.PFC_muFR_CFD_C1, C2.PFC_muFR_CFD_C1 );

MDcmi_block1 = CMI( C1_MD.MD_muFR_CFD_C1, C2_MD.MD_muFR_CFD_C1 );
MDcmi_block2 = CMI( C1rep_MD.MD_muFR_CFD_C1, C2_MD.MD_muFR_CFD_C1 );


indC1 = cmi_block1 > 0;
indC2 = cmi_block1 < 0;

%%
figure(1); set(gcf,'color','w');
plot( C1.RPFC, C2.RPFC,'o' ); hold on;
plot( linspace(0,5,2000), linspace(0,5,2000));
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')

figure(2); set(gcf,'color','w');
plot( C1.RPFC, C1rep.RPFC,'o' ); hold on;
plot( linspace(0,5,2000), linspace(0,5,2000));
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')


%%

figure(3); set(gcf,'color','w');
plot( C1_MD.RMD, C2_MD.RMD,'o' ); hold on;
plot( linspace(0,5,2000), linspace(0,5,2000));
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')

figure(4); set(gcf,'color','w');
plot( C1_MD.RMD, C1rep_MD.RMD,'o' ); hold on;
plot( linspace(0,5,2000), linspace(0,5,2000));
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')

%%

V = [PFC.TRRate];

VC2 = V(RangeC2, indC1); % context 1 preferring cells in C2

ind = Z_context2(:,10)==0;
Vnolas = VC2(ind, :);
Vlas   = VC2(~ind,:);

VC1r = V(RangeC1rep, indC1); % context 1 preferring cells in C1'

VC1 = V(RangeC1, indC1); % context 1 preferring cells in C1


%%

for i = 1:size(C1_MD.MD_dp_C1,1)
    CC(i) = corr(C1_MD.MD_dp_C1(i,:)', C1rep_MD.MD_dp_C1(i,:)');
    CC2(i) = corr(C1_MD.MD_dp_C1(i,:)', C2_MD.MD_dp_C1(i,:)');

end;

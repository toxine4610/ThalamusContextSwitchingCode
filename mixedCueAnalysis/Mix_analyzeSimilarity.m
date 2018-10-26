

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

[Spfc_mix_C1, Smd_mix_C1, ~, ~] = packageData(Z_context1, Smd, Spfc);
[Spfc_mix_C1rep, Smd_mix_C1rep, ~, ~] = packageData(Z_context1rep, Smd, Spfc);
[Spfc_mix_C2, Smd_mix_C2, ~, ~] = packageData(Z_context2, Smd, Spfc);



Out_C1 = extractVariables_PFC( Spfc_mix_C1 );
Out_C2 = extractVariables_PFC( Spfc_mix_C2 );
Out_C1re = extractVariables_PFC( Spfc_mix_C1rep );

%%

DP1 = Out_C1.PFC_dp_C1(goodPFC,:);
DP2 = Out_C1re.PFC_dp_C1(goodPFC,:);


for i = 1:size(DP1,1)
    foo1 = DP1(i,:);
    foo2 = DP2(i,:);
    CC(i) = corr( foo1', foo2' );
end
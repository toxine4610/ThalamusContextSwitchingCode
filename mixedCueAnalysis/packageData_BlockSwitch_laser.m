
function [Spfc_mix, Smd_mix, ZVis, ZAud, first, ErrorFrac] = packageData_BlockSwitch_laser(Z_C1, Smd, Spfc)


[~, ~, goodPFC, goodMD] = cleanData(Spfc, Smd, Z_C1);
Spfc =  Spfc(goodPFC == 1);
Smd  =  Smd(goodMD == 1);


%%

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

indCorrVis = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 1);
indCorrAud = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 2);

ZVis = Z_C1(indCorrVis, :);
ZAud = Z_C1(indCorrAud, :);

indRepvis = find( indCorrVis > SwitchTrial(2) );
indRepaud = find( indCorrAud > SwitchTrial(2) );

indFirstvis = find( indCorrVis <= SwitchTrial(1) );
indFirstaud = find( indCorrAud <= SwitchTrial(1) );

if Z_C1(2,9) == 0
    first = 0;
elseif Z_C1(2,9) == 3
    first = 3;
end;

indVisNoSound_NoLaser = find( ZVis(:,9) == 3 & ZVis(:,10) == 0);
indAudNoSound_NoLaser = find( ZAud(:,9) == 3 & ZAud(:,10) == 0);

indVisNoSound_Laser = find( ZVis(:,9) == 3 & ZVis(:,10) == 1);
indAudNoSound_Laser = find( ZAud(:, 9) == 3 & ZAud(:,10) == 1);

indVisSound_Laser = find( ZVis(:,9) == 0 & ZVis(:,10) == 1);
indAudSound_Laser = find( ZAud(:, 9) == 0 & ZAud(:,10) == 1);

indVisSound_NoLaser = find( ZVis(:,9) == 0 & ZVis(:,10) == 0);
indAudSound_NoLaser = find( ZAud(:, 9) == 0 & ZAud(:,10) == 0);

%
if first == 0
   indVisSound_NoLaser_first = intersect( indVisSound_NoLaser, indFirstvis );
   indAudSound_NoLaser_first = intersect( indAudSound_NoLaser, indFirstaud );
   
   indVisSound_NoLaser_rep = intersect( indVisSound_NoLaser, indRepvis );
   indAudSound_NoLaser_rep = intersect( indAudSound_NoLaser, indRepaud );
elseif first == 3
   indVisNoSound_NoLaser_first = intersect( indVisNoSound_NoLaser, indFirstvis );
   indAudNoSound_NoLaser_first = intersect( indAudNoSound_NoLaser, indFirstaud );
   
   indVisNoSound_NoLaser_rep = intersect( indVisNoSound_NoLaser, indRepvis );
   indAudNoSound_NoLaser_rep = intersect( indAudNoSound_NoLaser, indRepaud );
end;

indBadVis = find(Z_C1(:,1) == 0 & Z_C1(:,2) == 1);
indBadAud = find(Z_C1(:,1) == 0 & Z_C1(:,2) == 2);


ZVisBad = Z_C1(indBadVis, :);
ZAudBad = Z_C1(indBadAud, :);

indVisNoSoundBad_NoLaser = find( ZVisBad(:,9) == 3 & ZVisBad(:,10) == 0);
indAudNoSoundBad_NoLaser = find( ZAudBad(:, 9) == 3 & ZAudBad(:,10) == 0);

indVisNoSoundBad_Laser = find( ZVisBad(:,9) == 3 & ZVisBad(:,10) == 1);
indAudNoSoundBad_Laser = find( ZAudBad(:, 9) == 3 & ZAudBad(:,10) == 1);

indVisSoundBad_Laser = find( ZVisBad(:,9) == 0 & ZVisBad(:,10)== 1);
indAudSoundBad_Laser = find( ZAudBad(:, 9) == 0 & ZAudBad(:,10) == 1);

indVisSoundBad_NoLaser = find( ZVisBad(:,9) == 0 & ZVisBad(:,10)== 0);
indAudSoundBad_NoLaser = find( ZAudBad(:, 9) == 0 & ZAudBad(:,10) == 0);

%%

for i = 1:numel(Smd)
    
    if first == 0
    
        Smd_mix(i).SpikeTimes_R1_sound_NoLaser_first = Smd(i).SpikeTimes_R1C1(indVisSound_NoLaser_first);
        Smd_mix(i).SpikeTimes_R2_sound_NoLaser_first = Smd(i).SpikeTimes_R2C1(indAudSound_NoLaser_first);
    
        Smd_mix(i).SpikeTimes_R1_sound_NoLaser_rep = Smd(i).SpikeTimes_R1C1(indVisSound_NoLaser_rep);
        Smd_mix(i).SpikeTimes_R2_sound_NoLaser_rep = Smd(i).SpikeTimes_R2C1(indAudSound_NoLaser_rep);
    
        Smd_mix(i).SpikeTimes_R1_nosound_NoLaser = Smd(i).SpikeTimes_R1C1(indVisNoSound_NoLaser);
        Smd_mix(i).SpikeTimes_R2_nosound_NoLaser = Smd(i).SpikeTimes_R2C1(indAudNoSound_NoLaser);
    
        Smd_mix(i).SpikeTimes_R1_nosound_Laser = Smd(i).SpikeTimes_R1C1(indVisNoSound_Laser);
        Smd_mix(i).SpikeTimes_R2_nosound_Laser = Smd(i).SpikeTimes_R2C1(indAudNoSound_Laser);
    
    elseif first == 3
        
        Smd_mix(i).SpikeTimes_R1_sound_NoLaser = Smd(i).SpikeTimes_R1C1(indVisSound_NoLaser);
        Smd_mix(i).SpikeTimes_R2_sound_NoLaser = Smd(i).SpikeTimes_R2C1(indAudSound_NoLaser);
    
        Smd_mix(i).SpikeTimes_R1_sound_Laser = Smd(i).SpikeTimes_R1C1(indVisSound_Laser);
        Smd_mix(i).SpikeTimes_R2_sound_Laser = Smd(i).SpikeTimes_R2C1(indAudSound_Laser);
        
        Smd_mix(i).SpikeTimes_R1_nosound_NoLaser_first = Smd(i).SpikeTimes_R1C1(indVisNoSound_NoLaser_first);
        Smd_mix(i).SpikeTimes_R2_nosound_NoLaser_first = Smd(i).SpikeTimes_R2C1(indAudNoSound_NoLaser_first);
    
        Smd_mix(i).SpikeTimes_R1_nosound_NoLaser_rep = Smd(i).SpikeTimes_R1C1(indVisNoSound_NoLaser_rep);
        Smd_mix(i).SpikeTimes_R2_nosound_NoLaser_rep = Smd(i).SpikeTimes_R2C1(indAudNoSound_NoLaser_rep);
    end
end;

%%
for i = 1:numel(Spfc)
    
    if first == 0
    
        Spfc_mix(i).SpikeTimes_R1_sound_NoLaser_first = Spfc(i).SpikeTimes_R1C1(indVisSound_NoLaser_first);
        Spfc_mix(i).SpikeTimes_R2_sound_NoLaser_first = Spfc(i).SpikeTimes_R2C1(indAudSound_NoLaser_first);
    
        Spfc_mix(i).SpikeTimes_R1_sound_NoLaser_rep = Spfc(i).SpikeTimes_R1C1(indVisSound_NoLaser_rep);
        Spfc_mix(i).SpikeTimes_R2_sound_NoLaser_rep = Spfc(i).SpikeTimes_R2C1(indAudSound_NoLaser_rep);
    
        Spfc_mix(i).SpikeTimes_R1_nosound_NoLaser = Spfc(i).SpikeTimes_R1C1(indVisNoSound_NoLaser);
        Spfc_mix(i).SpikeTimes_R2_nosound_NoLaser = Spfc(i).SpikeTimes_R2C1(indAudNoSound_NoLaser);
    
        Spfc_mix(i).SpikeTimes_R1_nosound_Laser = Spfc(i).SpikeTimes_R1C1(indVisNoSound_Laser);
        Spfc_mix(i).SpikeTimes_R2_nosound_Laser = Spfc(i).SpikeTimes_R2C1(indAudNoSound_Laser);
    
    elseif first == 3
        
        Spfc_mix(i).SpikeTimes_R1_sound_NoLaser = Spfc(i).SpikeTimes_R1C1(indVisSound_NoLaser);
        Spfc_mix(i).SpikeTimes_R2_sound_NoLaser = Spfc(i).SpikeTimes_R2C1(indAudSound_NoLaser);
    
        Spfc_mix(i).SpikeTimes_R1_sound_Laser = Spfc(i).SpikeTimes_R1C1(indVisSound_Laser);
        Spfc_mix(i).SpikeTimes_R2_sound_Laser = Spfc(i).SpikeTimes_R2C1(indAudSound_Laser);
        
        Spfc_mix(i).SpikeTimes_R1_nosound_NoLaser_first = Spfc(i).SpikeTimes_R1C1(indVisNoSound_NoLaser_first);
        Spfc_mix(i).SpikeTimes_R2_nosound_NoLaser_first = Spfc(i).SpikeTimes_R2C1(indAudNoSound_NoLaser_first);
    
        Spfc_mix(i).SpikeTimes_R1_nosound_NoLaser_rep = Spfc(i).SpikeTimes_R1C1(indVisNoSound_NoLaser_rep);
        Spfc_mix(i).SpikeTimes_R2_nosound_NoLaser_rep = Spfc(i).SpikeTimes_R2C1(indAudNoSound_NoLaser_rep);
    end
end;

%%

indAudBlockLaser = find(Z_C1(:,10)==1 & Z_C1(:,9)==0);
indVisBlockLaser = find(Z_C1(:,10)==1 & Z_C1(:,9)==3);

indAudBlockNoLaser = find(Z_C1(:,10)==0 & Z_C1(:,9)==0);
indVisBlockNoLaser = find(Z_C1(:,10)==0 & Z_C1(:,9)==3);

foo  = Z_C1(indAudBlockLaser,:);
foo2 = Z_C1(indAudBlockNoLaser,:);
foo3 = Z_C1(indVisBlockNoLaser,:);
foo4  = Z_C1(indVisBlockLaser,:);

ErrorFrac(1) = length( find(foo(:,1)==0) )./size(foo,1);
ErrorFrac(2) = length( find(foo2(:,1)==0) )./size(foo2,1);
ErrorFrac(3) = length( find(foo3(:,1)==0) )./size(foo3,1);
ErrorFrac(4) = length( find(foo4(:,1)==0) )./size(foo4,1);

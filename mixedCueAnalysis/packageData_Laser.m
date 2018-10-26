

function [Spfc_mix, Smd_mix, ZVis, ZAud, first, ErrorFrac] = packageData_Laser(Z_C1, Smd, Spfc)


[~, ~, goodPFC, goodMD] = cleanData(Spfc, Smd, Z_C1);
Spfc =  Spfc(goodPFC == 1);
Smd  =  Smd(goodMD == 1);




indCorrVis = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 1);
indCorrAud = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 2);

ZVis = Z_C1(indCorrVis, :);
ZAud = Z_C1(indCorrAud, :);

if Z_C1(2,9) == 0
    first = 0;
elseif Z_C1(2,9) == 3
    first = 3;
end;

indVisNoSound_NoLaser = find( ZVis(:,9) == 3 & ZVis(:,10) == 0);
indAudNoSound_NoLaser = find( ZAud(:, 9) == 3 & ZAud(:,10) == 0);

indVisNoSound_Laser = find( ZVis(:,9) == 3 & ZVis(:,10) == 1);
indAudNoSound_Laser = find( ZAud(:, 9) == 3 & ZAud(:,10) == 1);

indVisSound_Laser = find( ZVis(:,9) == 0 & ZVis(:,10) == 1);
indAudSound_Laser = find( ZAud(:, 9) == 0 & ZAud(:,10) == 1);

indVisSound_NoLaser = find( ZVis(:,9) == 0 & ZVis(:,10) == 0);
indAudSound_NoLaser = find( ZAud(:, 9) == 0 & ZAud(:,10) == 0);

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
    
    Smd_mix(i).SpikeTimes_R1_sound_NoLaser = Smd(i).SpikeTimes_R1C1(indVisSound_NoLaser);
    Smd_mix(i).SpikeTimes_R2_sound_NoLaser = Smd(i).SpikeTimes_R2C1(indAudSound_NoLaser);
    
    Smd_mix(i).SpikeTimes_R1_nosound_NoLaser = Smd(i).SpikeTimes_R1C1(indVisNoSound_NoLaser);
    Smd_mix(i).SpikeTimes_R2_nosound_NoLaser = Smd(i).SpikeTimes_R2C1(indAudNoSound_NoLaser);
    
    Smd_mix(i).SpikeTimes_R1_sound_Laser = Smd(i).SpikeTimes_R1C1(indVisSound_Laser);
    Smd_mix(i).SpikeTimes_R2_sound_Laser = Smd(i).SpikeTimes_R2C1(indAudSound_Laser);
    
    Smd_mix(i).SpikeTimes_R1_nosound_Laser = Smd(i).SpikeTimes_R1C1(indVisNoSound_Laser);
    Smd_mix(i).SpikeTimes_R2_nosound_Laser = Smd(i).SpikeTimes_R2C1(indAudNoSound_Laser);
    
    
    Smd_mix(i).SpikeTimes_R1_sound_inc_NoLaser = Smd(i).SpikeTimes_R1C1_inc(indVisSoundBad_NoLaser);
    Smd_mix(i).SpikeTimes_R2_sound_inc_NoLaser = Smd(i).SpikeTimes_R2C1_inc(indAudSoundBad_NoLaser);
    
    Smd_mix(i).SpikeTimes_R1_nosound_inc_NoLaser = Smd(i).SpikeTimes_R1C1_inc(indVisNoSoundBad_NoLaser);
    Smd_mix(i).SpikeTimes_R2_nosound_inc_NoLaser = Smd(i).SpikeTimes_R2C1_inc(indAudNoSoundBad_NoLaser);
    
    Smd_mix(i).SpikeTimes_R1_sound_inc_Laser = Smd(i).SpikeTimes_R1C1_inc(indVisSoundBad_Laser);
    Smd_mix(i).SpikeTimes_R2_sound_inc_Laser = Smd(i).SpikeTimes_R2C1_inc(indAudSoundBad_Laser);
    
    Smd_mix(i).SpikeTimes_R1_nosound_inc_Laser = Smd(i).SpikeTimes_R1C1_inc(indVisNoSoundBad_Laser);
    Smd_mix(i).SpikeTimes_R2_nosound_inc_Laser = Smd(i).SpikeTimes_R2C1_inc(indAudNoSoundBad_Laser);
   
end;

%%
for i = 1:numel(Spfc)
    
    Spfc_mix(i).SpikeTimes_R1_sound_NoLaser = Spfc(i).SpikeTimes_R1C1(indVisSound_NoLaser);
    Spfc_mix(i).SpikeTimes_R2_sound_NoLaser = Spfc(i).SpikeTimes_R2C1(indAudSound_NoLaser);
    
    Spfc_mix(i).SpikeTimes_R1_nosound_NoLaser = Spfc(i).SpikeTimes_R1C1(indVisNoSound_NoLaser);
    Spfc_mix(i).SpikeTimes_R2_nosound_NoLaser = Spfc(i).SpikeTimes_R2C1(indAudNoSound_NoLaser);
    
    Spfc_mix(i).SpikeTimes_R1_sound_Laser = Spfc(i).SpikeTimes_R1C1(indVisSound_Laser);
    Spfc_mix(i).SpikeTimes_R2_sound_Laser = Spfc(i).SpikeTimes_R2C1(indAudSound_Laser);
    
    Spfc_mix(i).SpikeTimes_R1_nosound_Laser = Spfc(i).SpikeTimes_R1C1(indVisNoSound_Laser);
    Spfc_mix(i).SpikeTimes_R2_nosound_Laser = Spfc(i).SpikeTimes_R2C1(indAudNoSound_Laser);
    
    
    Spfc_mix(i).SpikeTimes_R1_sound_inc_NoLaser = Spfc(i).SpikeTimes_R1C1_inc(indVisSoundBad_NoLaser);
    Spfc_mix(i).SpikeTimes_R2_sound_inc_NoLaser = Spfc(i).SpikeTimes_R2C1_inc(indAudSoundBad_NoLaser);
    
    Spfc_mix(i).SpikeTimes_R1_nosound_inc_NoLaser = Spfc(i).SpikeTimes_R1C1_inc(indVisNoSoundBad_NoLaser);
    Spfc_mix(i).SpikeTimes_R2_nosound_inc_NoLaser = Spfc(i).SpikeTimes_R2C1_inc(indAudNoSoundBad_NoLaser);
    
    Spfc_mix(i).SpikeTimes_R1_sound_inc_Laser = Spfc(i).SpikeTimes_R1C1_inc(indVisSoundBad_Laser);
    Spfc_mix(i).SpikeTimes_R2_sound_inc_Laser = Spfc(i).SpikeTimes_R2C1_inc(indAudSoundBad_Laser);
    
    Spfc_mix(i).SpikeTimes_R1_nosound_inc_Laser = Spfc(i).SpikeTimes_R1C1_inc(indVisNoSoundBad_Laser);
    Spfc_mix(i).SpikeTimes_R2_nosound_inc_Laser = Spfc(i).SpikeTimes_R2C1_inc(indAudNoSoundBad_Laser);
   
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

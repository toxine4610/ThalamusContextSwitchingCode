function [Spfc_mix, Smd_mix, ZVis, ZAud] = packageData(Z_C1, Smd, Spfc)

[~, ~, goodPFC, goodMD] = cleanData(Spfc, Smd, Z_C1);
Spfc =  Spfc(goodPFC == 1);
Smd  =  Smd(goodMD == 1);


if size(Z_C1,2) == 9
    Z_C1(:,10) = 0;
end;

indCorrVis = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 1);
indCorrAud = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 2);

ZVis = Z_C1(indCorrVis, :);
ZAud = Z_C1(indCorrAud, :);

indVisNoSound = find( ZVis(:,9) == 3 & ZVis(:,10) == 0);
indAudNoSound = find( ZAud(:, 9) == 3 & ZAud(:,10) == 0);

indVisSound = find( ZVis(:,9) == 0 & ZVis(:,8)>=0.6 & ZVis(:,10) == 0);
indAudSound = find( ZAud(:, 9) == 0 & ZAud(:,8)>=0.6 & ZAud(:,10) == 0);

indBadVis = find(Z_C1(:,1) == 0 & Z_C1(:,2) == 1);
indBadAud = find(Z_C1(:,1) == 0 & Z_C1(:,2) == 2);

ZVisBad = Z_C1(indBadVis, :);
ZAudBad = Z_C1(indBadAud, :);

indVisNoSoundBad = find( ZVisBad(:,9) == 3& ZVisBad(:,10) == 0);
indAudNoSoundBad = find( ZAudBad(:, 9) == 3 & ZAudBad(:,10) == 0);

indVisSoundBad = find( ZVisBad(:,9) == 0& ZVisBad(:,10) == 0);
indAudSoundBad = find( ZAudBad(:, 9) == 0& ZAudBad(:,10) == 0);


%%

for i = 1:numel(Smd)
    
    Smd_mix(i).SpikeTimes_R1_sound = Smd(i).SpikeTimes_R1C1(indVisSound);
    Smd_mix(i).SpikeTimes_R2_sound = Smd(i).SpikeTimes_R2C1(indAudSound);
    
    Smd_mix(i).SpikeTimes_R1_nosound = Smd(i).SpikeTimes_R1C1(indVisNoSound);
    Smd_mix(i).SpikeTimes_R2_nosound = Smd(i).SpikeTimes_R2C1(indAudNoSound);
    
    Smd_mix(i).SpikeTimes_R1_sound_inc = Smd(i).SpikeTimes_R1C1_inc(indVisSoundBad);
    Smd_mix(i).SpikeTimes_R2_sound_inc = Smd(i).SpikeTimes_R2C1_inc(indAudSoundBad);
    
    Smd_mix(i).SpikeTimes_R1_nosound_inc = Smd(i).SpikeTimes_R1C1_inc(indVisNoSoundBad);
    Smd_mix(i).SpikeTimes_R2_nosound_inc = Smd(i).SpikeTimes_R2C1_inc(indAudNoSoundBad);
   
end;

for i = 1:numel(Spfc)
    Spfc_mix(i).SpikeTimes_R1_sound = Spfc(i).SpikeTimes_R1C1(indVisSound);
    Spfc_mix(i).SpikeTimes_R2_sound = Spfc(i).SpikeTimes_R2C1(indAudSound);
    
    Spfc_mix(i).SpikeTimes_R1_nosound = Spfc(i).SpikeTimes_R1C1(indVisNoSound);
    Spfc_mix(i).SpikeTimes_R2_nosound = Spfc(i).SpikeTimes_R2C1(indAudNoSound);
    
    Spfc_mix(i).SpikeTimes_R1_sound_inc = Spfc(i).SpikeTimes_R1C1_inc(indVisSoundBad);
    Spfc_mix(i).SpikeTimes_R2_sound_inc = Spfc(i).SpikeTimes_R2C1_inc(indAudSoundBad);
    
    Spfc_mix(i).SpikeTimes_R1_nosound_inc = Spfc(i).SpikeTimes_R1C1_inc(indVisNoSoundBad);
    Spfc_mix(i).SpikeTimes_R2_nosound_inc = Spfc(i).SpikeTimes_R2C1_inc(indAudNoSoundBad);

end;

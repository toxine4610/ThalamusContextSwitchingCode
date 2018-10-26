
indCorrVis = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 1 & Z_C1(:,10)==1);
indCorrAud = find(Z_C1(:,1) == 1 & Z_C1(:,2) == 2 & Z_C1(:,10)==1);

ZVis = Z_C1(indCorrVis, :);
ZAud = Z_C1(indCorrAud, :);

indVisNoSound_RightChoice = find( ZVis(:,9) == 3 & ZVis(:,5)==1);
indVisNoSound_LeftChoice = find( ZVis(:,9) == 3 & ZVis(:,5)==0);

indAudNoSound_RightChoice = find( ZAud(:, 9) == 3 & ZAud(:,5)==1);
indAudNoSound_LeftChoice = find( ZAud(:, 9) == 3 & ZAud(:,5)==0);

indVisSound_RightChoice  = find( ZVis(:,9) == 0 & ZVis(:,5)== 1);
indVisSound_LeftChoice  = find( ZVis(:,9) == 0 & ZVis(:,5)== 0);

indAudSound_RightChoice   = find( ZAud(:, 9) == 0 & ZAud(:,5)== 1);
indAudSound_LeftChoice   = find( ZAud(:, 9) == 0 & ZAud(:,5)== 0);

%%

range = [ -0.1, 1.2 ];
bin = 0.0010;
filtWidth = 0.08;


for i = 1:numel(Spfc)
    
    Spfc_mix(i).SpikeTimes_R1_sound_R = Spfc(i).SpikeTimes_R1C1(indVisSound_RightChoice);
    [R1_sound_R(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R1_sound_R, range , bin, filtWidth);

    Spfc_mix(i).SpikeTimes_R1_sound_L = Spfc(i).SpikeTimes_R1C1(indVisSound_LeftChoice);
    [R1_sound_L(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R1_sound_L, range , bin, filtWidth);

    Spfc_mix(i).SpikeTimes_R2_sound_R = Spfc(i).SpikeTimes_R2C1(indAudSound_RightChoice);
    [R2_sound_R(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R2_sound_R, range , bin, filtWidth);

    Spfc_mix(i).SpikeTimes_R2_sound_L = Spfc(i).SpikeTimes_R2C1(indAudSound_LeftChoice);
    [R2_sound_L(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R2_sound_L, range , bin, filtWidth);

    Spfc_mix(i).SpikeTimes_R1_nosound_R = Spfc(i).SpikeTimes_R1C1(indVisNoSound_RightChoice);
    [R1_nosound_R(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R1_nosound_R, range , bin, filtWidth);

    Spfc_mix(i).SpikeTimes_R1_nosound_L = Spfc(i).SpikeTimes_R1C1(indVisNoSound_LeftChoice);
    [R1_nosound_L(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R1_nosound_L, range , bin, filtWidth);
   
    Spfc_mix(i).SpikeTimes_R2_nosound_R = Spfc(i).SpikeTimes_R2C1(indAudNoSound_RightChoice);
    [R2_nosound_R(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R2_nosound_R, range , bin, filtWidth);

    Spfc_mix(i).SpikeTimes_R2_nosound_L = Spfc(i).SpikeTimes_R2C1(indAudNoSound_LeftChoice);
    [R2_nosound_L(i,:), ~, time, ~] = makeSpikeRates(Spfc_mix(i).SpikeTimes_R2_nosound_L, range , bin, filtWidth);

end;

%%
firingRatesAverage(:, 1, 1, :) = R1_nosound_L;
firingRatesAverage(:, 1, 2, :) = R1_nosound_R;

firingRatesAverage(:, 2, 1, :) = R2_nosound_L;
firingRatesAverage(:, 2, 2, :) = R2_nosound_R;

firingRatesAverage(:, 3, 1, :) = R1_sound_L;
firingRatesAverage(:, 3, 2, :) = R1_sound_R;

firingRatesAverage(:, 4, 1, :) = R2_sound_L;
firingRatesAverage(:, 4, 2, :) = R2_sound_R;

%%

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

timeEvents = [0, 1];

function  keep = cellCleaner(foo);

[FF_1, mean_1, var_1, bins, SC_1] = computeFanoFactor(foo.SpikeTimes_R1_nosound);
fract1 = findFractionBlankTrials(SC_1, bins);

[FF_2, mean_2, var_1, bins, SC_2] = computeFanoFactor(foo.SpikeTimes_R2_nosound);
fract2 = findFractionBlankTrials(SC_1, bins);

[FF_1, mean_1, var_1, bins, SC_3] = computeFanoFactor(foo.SpikeTimes_R1_sound);
fract3 = findFractionBlankTrials(SC_3, bins);

[FF_2, mean_2, var_1, bins, SC_4] = computeFanoFactor(foo.SpikeTimes_R2_sound);
fract4 = findFractionBlankTrials(SC_4, bins);

if (fract1 <= 0.25 && fract2 <= 0.25) || (fract3 <= 0.25 && fract4 <= 0.25)
    keep = 1;
else 
    keep = 0;
end;
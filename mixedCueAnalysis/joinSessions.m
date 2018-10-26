
A = load('C:\Users\Halassalab-CG\Dropbox\Rajeev\ForContextSwitchProject\DoubleCueDatabase\SOMCre\DualModalityCue\2017-12-20\Somcre_mixed_part1.mat');
B = load('C:\Users\Halassalab-CG\Dropbox\Rajeev\ForContextSwitchProject\DoubleCueDatabase\SOMCre\DualModalityCue\2017-12-20\Somcre_mixed_part2.mat');


Z_C1 = [A.Z_C1;B.Z_C1];

for i = 1:numel(A.Smd)
    Smd(i).filename = A.Smd(i).filename;
    Smd(i).unitNbr  = A.Smd(i).unitNbr;
    Smd(i).ttNbr    = A.Smd(i).ttNbr;
    Smd(i).SpikeTimes_R1C1 = [A.Smd(i).SpikeTimes_R1C1, B.Smd(i).SpikeTimes_R1C1];
    Smd(i).SpikeTimes_R2C1 = [A.Smd(i).SpikeTimes_R2C1, B.Smd(i).SpikeTimes_R2C1];
    Smd(i).SpikeTimes_R1C1_inc = [A.Smd(i).SpikeTimes_R1C1_inc, B.Smd(i).SpikeTimes_R1C1_inc];
    Smd(i).SpikeTimes_R2C1_inc = [A.Smd(i).SpikeTimes_R2C1_inc, B.Smd(i).SpikeTimes_R2C1_inc];
end;

for i = 1:numel(A.Spfc)
    Spfc(i).filename = A.Spfc(i).filename;
    Spfc(i).unitNbr  = A.Spfc(i).unitNbr;
    Spfc(i).ttNbr    = A.Spfc(i).ttNbr;
    Spfc(i).SpikeTimes_R1C1 = [A.Spfc(i).SpikeTimes_R1C1, B.Spfc(i).SpikeTimes_R1C1];
    Spfc(i).SpikeTimes_R2C1 = [A.Spfc(i).SpikeTimes_R2C1, B.Spfc(i).SpikeTimes_R2C1];
    Spfc(i).SpikeTimes_R1C1_inc = [A.Spfc(i).SpikeTimes_R1C1_inc, B.Spfc(i).SpikeTimes_R1C1_inc];
    Spfc(i).SpikeTimes_R2C1_inc = [A.Spfc(i).SpikeTimes_R2C1_inc, B.Spfc(i).SpikeTimes_R2C1_inc];
end;

D  = B.D;
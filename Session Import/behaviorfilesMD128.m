function [behaviorfilename]=behaviorfilesMD128(sn)
    BFN{1}.C1='Z:\CheetahData\MD128\2017-04-21_11-12-49\Behavior\MD128_neuralynx_042117.txt';
    BFN{1}.C2='Z:\CheetahData\MD128\2017-04-21_11-12-49\Behavior\MD128_neuralynx3Poke_042117.txt';
    BFN{2}.C1='Z:\CheetahData\MD128\2017-04-25_16-57-55\behaviour\MD128_neuralynx_042517.txt';
    BFN{2}.C2='Z:\CheetahData\MD128\2017-04-25_16-57-55\behaviour\MD128_neuralynx3Poke_042517.txt';
    BFN{3}.C1='Z:\CheetahData\MD128\2017-04-26_13-49-09\behaviour\MD128_neuralynx_042617.txt';
    BFN{3}.C2='Z:\CheetahData\MD128\2017-04-26_13-49-09\behaviour\MD128_neuralynx3Poke_042617.txt';
    BFN{4}.C1='';
    BFN{4}.C2='';
    BFN{5}.C1='';
    BFN{5}.C2='';
    BFN{6}.C1='';
    BFN{6}.C2='';
    BFN{7}.C1='Y:\CheetahData\MD128\2017-05-27_14-26-56\behavior\MD128_C1_3pokes_20170527.txt';
    BFN{7}.C2='Y:\CheetahData\MD128\2017-05-27_14-26-56\behavior\MD128_C2_2-4pokes_20170527.txt';
    behaviorfilename.C1=BFN{sn}.C1;
    behaviorfilename.C2=BFN{sn}.C2;
end
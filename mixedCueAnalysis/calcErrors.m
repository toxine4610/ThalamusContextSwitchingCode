
function [pCorr_given_2_C1, pWrong_given_2_C1, C_C1] = calcErrors(Z)

corr_given_2 = []; wrong_given_2 = [];


for i = 3:size(Z,1)
    
    curr    = Z(i,1);
    curr_1  = Z(i-1,1);
    curr_2  = Z(i-2,1);
    
    if curr_1 == 1 && curr_2 == 1 && curr == 1
       corr_given_2 = [corr_given_2, i];
    elseif  curr_1 == 1 && curr_2 == 1 && curr == 0
       wrong_given_2 = [wrong_given_2, i];
    end;
    
end;

pCorr_given_2_C1 = length(corr_given_2)./( length(corr_given_2) + length(wrong_given_2) );
pWrong_given_2_C1 = 1-pCorr_given_2_C1;
C_C1 = length(corr_given_2);
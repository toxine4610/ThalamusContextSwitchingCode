
Ztruncated = Z_C1( :, :);


indNoLas = find( Ztruncated(:,10) == 0 & Ztruncated(:,9) == 0);
Z        = Ztruncated(indNoLas,:);

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

indNoLas = find( Ztruncated(:,10) == 0 & Ztruncated(:,9) == 3);
Z        = Ztruncated(indNoLas,:);

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

pCorr_given_2_C2 = length(corr_given_2)./( length(corr_given_2) + length(wrong_given_2) );
pWrong_given_2_C2 = 1-pCorr_given_2_C2;
C_C2 = length(corr_given_2);



figure(100); set(gcf,'color','w');
subplot(1,2,1)
b=bar(1, pCorr_given_2_C1); hold on
b(1).FaceColor =  'k'; b(1).BarWidth = 0.4;
axis square; box off; set(gca,'tickdir','out','fontsize',16);
title('Context 1'); ylabel('Fraction Correct');

subplot(1,2,2)
b=bar(1, pCorr_given_2_C2); hold on
b(1).FaceColor =  'k';b(1).BarWidth = 0.4;
axis square; box off; set(gca,'tickdir','out','fontsize',16);
title('Context 2');ylabel('Fraction Correct');

%%
indLas   = find( Ztruncated(:,10) == 1 & Ztruncated(:,9) == 0);
Z        = Ztruncated(indLas,:);
numCorrec  = length( find( Z(:,1) == 1) );
pCorrect_C1 = numCorrec./size(Z,1);
C_C1_las = numCorrec;

indLas   = find( Ztruncated(:,10) == 1 & Ztruncated(:,9) == 3);
Z        = Ztruncated(indLas,:);
numCorrec  = length( find( Z(:,1) == 1) );
pCorrect_C2 = numCorrec./size(Z,1);
C_C2_las = numCorrec;


figure(100); set(gcf,'color','w');
subplot(1,2,1)
b=bar(1.4, pCorrect_C1); hold on
b(1).FaceColor =  'y';b(1).BarWidth = 0.4;
ylim([0,1]);

subplot(1,2,2);
b=bar(1.4, pCorrect_C2); hold on
b(1).FaceColor =  'y';b(1).BarWidth = 0.4;
axis square; box off; set(gca,'tickdir','out','fontsize',16);
ylim([0,1]); 

%%

x = table([C_C1;C_C2],[C_C1_las;C_C2_las],'VariableNames',{'NoLas','Las'},'RowNames',{'C1','C2'})
[h,p,stats] = fishertest(x,'Tail','both','Alpha',0.05)


% addpath C:\Users\Halassalab-CG\Dropbox (Personal)\Rajeev\GLM\mdpfcglm

dataRep = '../ForContextSwitchProject/DoubleCueDatabase/BPM4_3/BlockSwitch/';
Dates = getDateList(dataRep);

%%

for d = setdiff(4:12,[6,11])
    
    dataName = [ Dates{d} filesep 'BPM43_block_withlaser_uni_Session.mat'];
    load( [dataRep dataName] );
 
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
%     [pCorr_given_2_C1, pWrong_given_2_C1, C_C1] = calcErrors(Z_context1);
%     [pCorr_given_2_C1re, pWrong_given_2_C1re, C_C1re] = calcErrors(Z_context1rep);
%     
%     indlaser = Z_context2(:,10) == 1;
%     [pCorr_given_2_C2las, pWrong_given_2_C2, C_C2] = calcErrors(Z_context2(indlaser, :));
%     
%     indlaser = Z_context2(:,10) == 0;
%     [pCorr_given_2_C2nolas, pWrong_given_2_C2, C_C2] = calcErrors(Z_context2(indlaser, :));
%     
%     
%     figure(1);
%     
%     b=bar(1, pCorr_given_2_C1); hold on
%     b(1).FaceColor =  'k';
%     
%     b=bar(2, pCorr_given_2_C2nolas); hold on
%     b(1).FaceColor =  'r';
%     
%     % b=bar(3, pCorr_given_2_C2las); hold on
%     % b(1).FaceColor =  'y';
%     b=bar(3, pCorr_given_2_C1re); hold on
%     b(1).FaceColor =  'k';
%     
%     axis square; box off; set(gca,'tickdir','out','fontsize',16);
%     ylabel('Fraction Correct');
    
    
    %% plots black line -- real data
   
    w = gausswin(5);
    w = w/sum(w);
    
    perf1 =conv(Z_C1(:,1),w,'same');
    perf2 =conv(Z_C2(:,1),w,'same');
%     perf3 =conv(Z_context1rep(:,1),w,'same');
    
    P = [smooth(perf1) ; smooth(perf2)];
    
    figure(d); set(gcf,'color','w');
    plot( 1:length(P), P, 'color','k','linewidth',3 ); hold on
    for i = 1:length(SwitchTrial)
        line( [SwitchTrial(i), SwitchTrial(i)], [0,1],'color','r' );
    end;
    
    line([0,200],[0.5, 0.5],'color','r')
    axis normal; box off; set(gca,'tickdir','out','fontsize',16);
    
%     [c_data, t_data] = hpfilter( P, 1600 );
%     plot( 1:length(P), c_data );
    
    %% plots model -- blue line.
    
    [ExpH_C1, rup_C1, rlo_C1] = computeExpectation(Z_C1);
    [ExpH_C2, rup_C2, rlo_C2] = computeExpectation(Z_C2);
    [ExpH_C3, rup_C3, rlo_C3] = computeExpectation(Z_context1rep);
    [ExpH_All, rup, rlo] = computeExpectation(Z_C1);
    
    figure(d);
    
%     yyaxis left;
    plot( 1:length(ExpH_C1), (ExpH_C1), 'linewidth',3); hold on
    plot( linspace(SwitchTrial(1), SwitchTrial(2), length(ExpH_C2)), (ExpH_C2),'linewidth',3 );
    plot( linspace(SwitchTrial(2), max(size(Z_C1)), length(ExpH_C3)), (ExpH_C3),'linewidth',3 );
    plot( 1:size(Z_C1,1), smooth(ExpH_All),'r','linewidth',3);
    
%     yyaxis right;
%     stem( 1:size(Z_C1,1), 2*Z_C1(:,1)-1 ); 
    
    
%     CC_whole(d-2) = sum( (ExpH_All'-mean(P)).^2 )./ sum( (P-mean(P)).^2 );
%     CC_B1(d-2) = sum( (ExpH_C1'-mean(perf1)).^2 )./ sum( (perf1-mean(perf1)).^2 );
%     CC_B2(d-2) = sum( (ExpH_C2'-mean(perf2)).^2 )./ sum( (perf2-mean(perf2)).^2 );
%     CC_B3(d-2) = sum( (ExpH_C3'-mean(perf3)).^2 )./ sum( (perf3-mean(perf3)).^2 );
    
    CC_whole(d-2) = corr( ExpH_All', P );
    CC_B1(d-2)    = corr( ExpH_C1', (perf1));
    CC_B2(d-2)    = corr( ExpH_C2', (perf2));
    CC_B3(d-2)    = corr( ExpH_C3', (perf3));
% %     %%
%     Tc = [];
%     y  = smooth(perf1);
%     for i = 1:length( y )
%         if y(i) >= ExpH_C1(i)
%             Yconsistent(i) = y(i);
%             Tc = [Tc,i];
%         else
%             Yconsistent(i) = NaN;
%         end;
%     end;
%     
%     Z_context1_consistent = Z_context1( Tc, : );
%     [pCorr_given_2_C1_corrected, ~, ~] = calcErrors(Z_context1_consistent);
%     
%     
%     Tc = [];
%     y  = smooth(perf2);
%     for i = 1:length( y )
%         if y(i) >= ExpH_C2(i)
%             Yconsistent(i) = y(i);
%             Tc = [Tc,i];
%         else
%             Yconsistent(i) = NaN;
%         end;
%     end;
%     
%     Z_context2_consistent = Z_context2( Tc, : );
%     [pCorr_given_2_C2_corrected, ~, ~] = calcErrors(Z_context2_consistent);
%     C1_C2_Transition = Tc;
%     
%     
%     Tc = [];
%     y  = smooth(perf3);
%     for i = 1:length( y )
%         if y(i) >= ExpH_C3(i)
%             Yconsistent(i) = y(i);
%             Tc = [Tc,i];
%         else
%             Yconsistent(i) = NaN;
%         end;
%     end;
%     
%     Z_context1rep_consistent = Z_context1rep( Tc, : );
%     C2_C1_Transition = Tc;
%     [pCorr_given_2_C1re_corrected, ~, ~] = calcErrors(Z_context2_consistent);
%     
%     
%     
%     figure(3);
%     
%     b=bar(1, pCorr_given_2_C1_corrected); hold on
%     b(1).FaceColor =  'k';
%     
%     b=bar(2, pCorr_given_2_C2_corrected); hold on
%     b(1).FaceColor =  'r';
%     
%     b=bar(3, pCorr_given_2_C1re_corrected); hold on
%     b(1).FaceColor =  'k';
%     
%     plot(1, mean(ExpH_C1),'o','markersize',12,'markerfacecolor','g');
%     plot(2, mean(ExpH_C2),'o','markersize',12,'markerfacecolor','g');
%     plot(3, mean(ExpH_C3),'o','markersize',12,'markerfacecolor','g');
%     
%     plot(1, mean(rlo_C1(2:end)),'o','markersize',12,'markerfacecolor','g');
%     plot(2, mean(rlo_C2(2:end)),'o','markersize',12,'markerfacecolor','g');
%     plot(3, mean(rlo_C3(2:end)),'o','markersize',12,'markerfacecolor','g');
%     axis square; box off; set(gca,'tickdir','out','fontsize',16);
%     ylabel('Fraction Correct');
%     ylim([0,1])


end;



%%

cmap = lines(3);

for i = 1:length(CC_whole)
    figure(400); set(gcf,'color','w')
    plot( i, CC_whole(i),'o','markerfacecolor','r','markersize',15); hold on
    plot( i, CC_B1(i),'o', 'markerfacecolor', cmap(1,:),'markersize',15);
    plot( i, CC_B2(i),'o', 'markerfacecolor', cmap(2,:),'markersize',15);
    plot( i, CC_B3(i),'o', 'markerfacecolor', cmap(3,:),'markersize',15);
end
plot( 1:length(CC_whole), CC_whole,'r')
plot( 1:length(CC_whole), CC_B1,'color',cmap(1,:))
plot( 1:length(CC_whole), CC_B2,'color',cmap(2,:))
plot( 1:length(CC_whole), CC_B3,'color',cmap(3,:))
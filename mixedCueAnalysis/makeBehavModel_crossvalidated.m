


numCVSplits  = 5;
history = 5;
Indices = crossvalind( 'Kfold', size(Vmega,1) , numCVSplits) ;


%%
for i = 1: numCVSplits
    
    clear Vtrain Vtest X Y X_TEST Y_TEST Xnew Xnew_TEST B
    
    
    Vtrain = Vmega( Indices ~= i, :);
    Vtest  = Vmega( Indices == i, :);
    
    %% ===== make regressors...............................................
    [currRule, currRew, currChoice, currContext, currLaser, R, RW, N, NW, L, C ] = getRegressors(Vtrain, history);
    Y = categorical(currChoice);
    X = [currRule', currRew', (currLaser.*currRule)', R, N, RW, NW, L.*R, L.*N, L.*RW, L.*NW] ;
    
%     [Xnew, ind] = cleanDesignMatrix(X);
    ind = []
    [B, Borig ,dev, stats] = performRegression(X, Y, ind, 1 );
    
    CV(i).Betanew = B;
    CV(i).Beta = Borig;
    CV(i).Deviance = dev;
    CV(i).Vtest  = Vtest;
    CV(i).Vtrain = Vtrain;
    
    %% ==== evaluate regression............................................
    clear X_TEST Xnew_TEST
    clear currRule currRew currChoice currContext currLaser R RW N NW L C
    [currRule, currRew, currChoice, currContext, currLaser, R, RW, N, NW, L, C ] = getRegressors(Vtest, history);
    Y_TEST = categorical(currChoice);
    X_TEST = [currRule', currRew', (currLaser.*currRule)', R, N, RW, NW, L.*R, L.*N, L.*RW, L.*NW] ;
    
%     Xnew_TEST = X_TEST( :, setdiff(1:size(X_TEST,2), ind));
% 
%     
    [pihat,dlow,dhi] = mnrval(B,X_TEST,stats);
    
    [r(i), p(i)] = corr( currChoice', pihat(:,2),'type','Kendall');
    
    figure(i)
    yyaxis right
    plot( 1:length(currChoice), smooth( pihat(:,2), 2 ) ); hold on;
    yyaxis left
    plot( 1:length(currChoice), currChoice ); ylim([-2,2]); hold on
end

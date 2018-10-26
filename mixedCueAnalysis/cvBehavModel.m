



for iter = 1:200
    history = 10;
    
    numCVSplits  = 5;
    Indices = crossvalind( 'Kfold', size(Vmega,1) , numCVSplits) ;
    
    for cv = 1:numCVSplits
        
        VmegaCV = [];
        
        test = (Indices == cv);
        train = ~test;
        
        % ========= train model
        VmegaCV_train = Vmega( train,: );
        
        % clear B dev stats Y X
        [currRule, currRew, currChoice, currContext, currLaser, ....
            R, RW, N, NW, L, C ] = getRegressors(VmegaCV_train, history);
        Y = categorical(currChoice);
        X = [currRule', currRew', (currLaser.*currRule)', R, N, RW, NW, L.*R, L.*N, L.*RW, L.*NW] ;
        [B, dev, stats] = mnrfit( X, Y,'Interactions','on' );
        
        % ======== test model
        
        clear X Y
        
        VmegaCV_test = Vmega( test,: );
        
        [currRule, currRew, currChoice, currContext, currLaser, ....
            R, RW, N, NW, L, C ] = getRegressors(VmegaCV_test, history);
        Y = categorical(currChoice);
        X = [currRule', currRew', (currLaser.*currRule)', R, N, RW, NW, L.*R, L.*N, L.*RW, L.*NW] ;
        
        [pihat,dlow,dhi] = mnrval(B,X,stats);
        
        pihatdisc( pihat(:,2) > 0.5) = 1;
        pihatdisc( pihat(:,2) <= 0.5) = -1;
        
        r(iter,cv) =  corr( pihatdisc', currChoice'  );
        
        CV(cv).Beta = B;
        CV(cv).stats = stats;
    end;
    
end;
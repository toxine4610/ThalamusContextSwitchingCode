
PTRC2Vec = [];
PTRVec   = [];


for i = 1:9
    A = load( ['MX_Sess_' num2str(i) '_C1.mat'] );
    B = load( ['MX_Sess_' num2str(i) '_C2.mat'] );
    
    llc = A.Lc./1e4;
    indkeep = llc > 0.3;
    
    for i = 1:size(A.CouplingFilterPFC,2)
        foo = A.CouplingFilter{i};
        for j = 1:size(foo,2)
            properties(j) = extractProps( A.CouplingFilterPFC{i}{j} );
        end;
        PTR(i,:) = [properties.PTR];
        Rat(i,:) = [properties.Diff];
        Ratio(i,:) = [properties.Ratio];
    end;
    
    for i = 1:size(B.CouplingFilterPFC,2)
        foo = B.CouplingFilter{i};
        for j = 1:size(foo,2)
            properties(j) = extractProps( B.CouplingFilterPFC{i}{j} );
        end;
        PTRC2(i,:) = [properties.PTR];
        RatC2(i,:) = [properties.Diff];
        RatioC2(i,:) = [properties.Ratio];
    end;
    
    Rat = Rat(indkeep,:);
    RatC2 = RatC2(indkeep,:);
    
    PTR  =PTR(indkeep,:);
    PTRC2 = PTRC2(indkeep,:);
    
    PTRVec = [PTRVec; median(PTR,2)];
    PTRC2Vec = [PTRC2Vec;median(PTRC2,2)];
    
    %%
    figure(1);
    plot( median(Rat,2), median(RatC2,2),'o','markersize',10,'color','k'); hold on;
    plot( linspace(-3,1,2000), linspace(-3,1,2000) );
    
    
    figure(2);
    plot( median(PTRC2,2), median(PTR,2),'o','markersize',10,'color','k'); hold on;
    plot( linspace(0,3,2000), linspace(0,3,2000) );
end
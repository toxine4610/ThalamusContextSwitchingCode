X = cell2mat( AUC' );
X( X < -10 ) = NaN;
X( X > 10 )  = NaN;

poorFits =  find(Lc < 1e3 );
X(poorFits,:) = NaN;

SMD = sign( CMIMD );

for i = 1:numPFC
    
    cmi_this = CMIPFC(i);
    s = sign(cmi_this);
    
    ind = find( SMD == s);
    IndexSame{i} = ind;
    ind2 = find( SMD ~= s);
    IndexOpp{i}  = ind2;
    
    ct = 0;
    for j = ind
        ct = ct+1;
        dist{i}(ct) = sqrt( cmi_this.^2 + CMIMD(j).^2 );
    end;
    
    ct = 0;
    for j = ind2
        ct = ct+1;
        distOpp{i}(ct) = sqrt( cmi_this.^2 + CMIMD(j).^2 );
    end;
    
end;

%%

for i = 1:numPFC
    
    cmi_this = CMIPFC(i);
    
    for j = 1:numMD
        distAll(i,j) = sqrt( cmi_this.^2 + CMIMD(j).^2 );
    end;

end;

%%

CFS_same = [];
for i = 1:numPFC
    CFS_same = [CFS_same, X(i, IndexSame{i} )];
end;

CFS_oppo = [];
for i = 1:numPFC
    CFS_oppo = [CFS_oppo, X(i, IndexOpp{i} )];
end;

%%

for i = setdiff( 1:numPFC, poorFits );
    figure(1); set(gcf,'color','w');
    plot( dist{i}, X(i, IndexSame{i}),'o','markersize',12','markerfacecolor','k','markeredgecolor','w');
    hold on
end;
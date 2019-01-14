X = cell2mat( A.AUC' );
lls  = A.Lc./1e4;
indok = find(lls > 0.3);


CFC1 = A.CouplingFilter;
CFC2 = B.CouplingFilter;

Y = [];

for i = 1:length(indok)
    
    cmidist_thiscell = dist( indok(i), : );
    mdin_thiscell = X(indok(i),:);
    
    CF_thisCell =  CFC2{ indok(i) };
    for j = 1:size(CF_thisCell,2)
        propertiesC2(j) = extractProps(CF_thisCell{j});
    end;
    
    CF_thisCell =  CFC1{ indok(i) };
    for j = 1:size(CF_thisCell,2)
        propertiesC1(j) = extractProps(CF_thisCell{j});
    end;
    
    PTR = [propertiesC1.PTR];
    Ratio = [propertiesC1.Ratio];
    EField = [propertiesC1.ExcSubfield];
    IField = [propertiesC1.InhibSubfield];
    
    
    PTRC2 = [propertiesC2.PTR];
    RatioC2 = [propertiesC2.Ratio];
    EFieldC2 = [propertiesC2.ExcSubfield];
    IFieldC2 = [propertiesC2.InhibSubfield];
    
    
    figure(1);
    plot( mean( cmidist_thiscell(cmidist_thiscell<0) ), mean( PTR(cmidist_thiscell<0) ),'o','color','r'); hold on;
    plot( mean( cmidist_thiscell(cmidist_thiscell>0) ), mean( PTR(cmidist_thiscell>0) ),'o','color','b'); hold on;
    
    
    figure(2);
    plot( median( cmidist_thiscell(cmidist_thiscell<0) ), median( Ratio(cmidist_thiscell<0) ),'o','color','r'); hold on;
    plot( median( cmidist_thiscell(cmidist_thiscell>0) ), median( Ratio(cmidist_thiscell>0) ),'o','color','b'); hold on;
    

    figure(3);
    plot( ( cmidist_thiscell(cmidist_thiscell<0) ), ( EField(cmidist_thiscell<0) ),'o','color','r'); hold on;
    plot( ( cmidist_thiscell(cmidist_thiscell>0) ), ( EField(cmidist_thiscell>0) ),'o','color','b'); hold on;
    
    figure(4);
    plot( median(PTR), median(PTRC2), 'o', 'color','k'); hold on
    ptr{i} = PTR;
    ptrc2{i} = PTRC2;
end;

    
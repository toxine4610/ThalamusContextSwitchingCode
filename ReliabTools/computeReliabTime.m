
function [ES,ED] = computeReliabTime( nBins, binEdges, windowLen, windowType, SpikeTimesCell, P)


%%
for b = 1 : nBins
    
    switch windowType
        case 'Sliding'
            bin1  = binEdges(b);
            bin2  = binEdges(b)+windowLen;
        case 'Fixed'
            bin1  = timeBins(b);
            bin2  = timeBins(b+1);
    end;
    
    for i = 1:size(P,1)
        
        c1 = P(i,1); c2 = P(i,2);
        
        sTrain1 = SpikeTimesCell{c1}( SpikeTimesCell{c1} >= bin1 &  SpikeTimesCell{c1} <= bin2  );
        sTrain2 = SpikeTimesCell{c2}( SpikeTimesCell{c2} >= bin1 &  SpikeTimesCell{c2} <= bin2  );
        
        try
            if sTrain2 < 1e-3
                sTrain2 = [];
            end;
            
            if sTrain1 < 1e-3
                sTrain1 = [];
            end;
            
            s = find(sTrain1 == 0);
            if ~isempty(s)
                sTrain1(s) = sTrain1(s)+0.001 ;
            end;
            
            s = find(sTrain2 == 0);
            if ~isempty(s)
                sTrain2(s) = sTrain2(s)+0.001 ;
            end;
            
            
            if ~isempty(sTrain1) && ~isempty(sTrain2)
                
                [es(i),ed(i)] = Event_Sync( sTrain1, sTrain2 );
            else
                es(i) = NaN;
                ed(i) = NaN;
            end;
        catch err
            es(i) = NaN;
            ed(i) = NaN;
        end;
    end;
    
    ES{b} = es( ~isnan(es) );
    ED{b} = ed( ~isnan(ed) );
    
    clear es ed
end;
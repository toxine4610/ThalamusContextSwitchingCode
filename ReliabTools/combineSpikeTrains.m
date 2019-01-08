function S = combineSpikeTrains(SpikeTimesCell_1T, SpikeTimesCell_2T);

sizeNew = max(size(SpikeTimesCell_1T)) + max(size(SpikeTimesCell_2T));
    
    for i = 1:sizeNew
        if i <= max(size(SpikeTimesCell_1T))
            S{i} = SpikeTimesCell_1T{i};
        else
            S{i} = SpikeTimesCell_2T{i - max( size(SpikeTimesCell_1T))};
        end;
    end;
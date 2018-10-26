
function props = computeSpikeStatistics(waveforms_this_unit);


% compute spike statistics ================================
[neg_pks, neg_loc] = findpeaks(-waveforms_this_unit,'MinPeakHeight',0,'SortStr','descend');
[pos_pks, pos_loc] = findpeaks(waveforms_this_unit,'MinPeakHeight',0,'SortStr','descend');


Peak = pos_pks(1);
Trough = neg_pks(1);

props.tpeak = pos_loc(1);
props.ttrough = neg_loc(1);
props.peaktroughratio = Peak./abs(Trough);
props.SpikeAmplitude   = Peak - (Trough);


% % % props.numNegPeaks = length(neg_pks.loc);
% % % props.numPosPeaks = length(pos_pks.loc);
% % % 
% % % for i = 1:length(neg_pks.loc)
% % %     NPk(i) = waveforms_this_unit( neg_pks.loc(i) );
% % % end;
% % % for i = 1:length(pos_pks.loc)
% % %     PPk(i) = waveforms_this_unit( pos_pks.loc(i) );
% % % end;
% % % 
% % % [Peak, tpeak] = max( PPk );
% % % props.tpeak         =  pos_pks.loc(tpeak);
% % % 
% % % [Trough, ttrough] = min( NPk );
% % % props.ttrough         =  neg_pks.loc(ttrough);




[pk, tmax] = max(waveforms_this_unit);
[tr, tmin] = min(waveforms_this_unit);


props.cycle = abs(tmax-tmin);

halfMax  = pk./2;
indhalfMax = find( waveforms_this_unit >= halfMax );
if ~isempty(indhalfMax)
    props.SpikeWidth = indhalfMax(end)-indhalfMax(1);
else
    props.SpikeWidth = NaN;
end;

props.PTTime    = abs(props.tpeak - props.ttrough);
% props.RepolTime = pos_pks.loc(end) - neg_pks.loc(end);
% props.RiseTime  = props.tpeak  - neg_pks.loc(1);

props.waveform  = waveforms_this_unit;



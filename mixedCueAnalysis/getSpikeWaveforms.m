

% function getSpikeWaveforms

%%

clear props;

ct = 0;

PFCTT = [1:12];

folder = uipickfiles('FilterSpec', 'X:\Rajeev\Cheetah\BPM1-3');
diagnostic = 1;


for i = 1:size(folder,2)
    %
    
    this_folder = folder{i};
    v = this_folder(11:end);
    fprintf('Processing: %s\n', this_folder);
    
    for n = 1:length(PFCTT)
        this_TT = ['TT' num2str(PFCTT(n)) 'cluster'];
        fprintf('\tLoading TT #%d...', n)
        A = load( [this_folder filesep this_TT '.mat'] );
        eval(['labels_thisTT = A.TT' num2str(PFCTT(n)) '.labels' ';']);
        eval(['waveforms_thisTT = A.TT' num2str(PFCTT(n)) '.waveforms' ';']);
        
        s{i}(n) = size(labels_thisTT,2);
        
        for u = 1:size(labels_thisTT,2)
            try
                ct = ct+1;
                labels_this_unit = labels_thisTT{u};
                waveforms_this_unit = squeeze( (waveforms_thisTT(:,:, labels_this_unit) ));
                
                waveform_this_unit  = svdDecompWaveforms(waveforms_this_unit);
                
                figure(100); plot(1:32, waveform_this_unit); hold on;
                axis square; box off; set(gca,'tickdir','out','fontsize',16);
                
                %                 waveforms_this_unit = smooth( median( waveforms_this_unit, 2) );
                %                 waveforms_this_unit  = interp1(1:32, waveforms_this_unit, linspace(0,32, 2000) );
                props(ct) = computeSpikeStatistics(waveform_this_unit);
                count(ct) = size(labels_this_unit,2);
            catch err
            end
        end;
        fprintf('...Done\n')
    end;
    
    
end;



%%
%

rng(1);


c = [props.cycle]
indx = find(c < 16);

PTR = [props.peaktroughratio];
PTT = [props.PTTime];
SW  = [props.SpikeWidth];
SA  = [props.SpikeAmplitude];

PTR = PTR(indx);
PTT = PTT(indx);
SW  = SW(indx);
SA  = SA(indx);

ind = find( PTR > prctile(PTR,2) & PTR < prctile(PTR,98) );

PTT = PTT(ind);
PTR = PTR(ind);
SW  = SW(ind);
SA  = SA(ind);

X = [ abs(SA); SW; PTT ; PTR];
opts = statset('Display','final');
[idx,C] = kmeans(X', 4,'Replicates', 25,'Distance','cityblock','Options',opts);

figure(2);
plot( PTT(idx==1), PTR(idx==1), 'o', 'markerfacecolor','r','markeredgecolor','w'); hold on;
plot( PTT(idx==2), PTR(idx==2), 'o', 'markerfacecolor','b','markeredgecolor','w'); hold on;
plot( PTT(idx==3), PTR(idx==3), 'o', 'markerfacecolor','g','markeredgecolor','w'); hold on;
plot( PTT(idx==4), PTR(idx==4), 'o', 'markerfacecolor','m','markeredgecolor','w'); hold on;

axis square; box off; set(gca,'tickdir','out','fontsize',16);


figure(3);
plot( PTT(idx==1), count(idx==1), 'o', 'markerfacecolor','r','markeredgecolor','w'); hold on;
plot( PTT(idx==2), count(idx==2), 'o', 'markerfacecolor','b','markeredgecolor','w'); hold on;
plot( PTT(idx==3), count(idx==3), 'o', 'markerfacecolor','g','markeredgecolor','w'); hold on;
plot( PTT(idx==4), count(idx==4), 'o', 'markerfacecolor','m','markeredgecolor','w'); hold on;

axis square; box off; set(gca,'tickdir','out','fontsize',16);

%%
W = {props.waveform};
W = cell2mat( W' );

% 
% figure(1); plot( (W(idx==1,:)' ),'r'); hold on
% figure(1); plot( (W(idx==2,:)' ),'b'); hold on
% figure(1); plot( (W(idx==3,:)' ),'g'); hold on
figure(1); plot( (W(idx==4,:)' ),'m'); hold on
axis square; box off; set(gca,'tickdir','out','fontsize',16);


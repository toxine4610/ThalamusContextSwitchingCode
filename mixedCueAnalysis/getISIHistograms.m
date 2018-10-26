

ct = 0;

MDTT = 10:16;

folder = uipickfiles('FilterSpec', 'F:\somcre');



%%
for i = 1:size(folder,2)
    
    try
        this_folder = folder{i};
        v = this_folder(11:end);
        fprintf('Processing: %s\n', this_folder);
        
        for T = MDTT
            
            fprintf('\tLoading TT #%d...', T)
            C1 = load([ this_folder '\GlobalpeakC1\' v '_PFCMDSection1TT' num2str(T) '.mat'] );
            fprintf('..Done\n');
            
            numunits = length( C1.num_seq );
            
            for u = 1:numunits
                clear spiketimes_thiscell N edges
                ct = ct+1;
                eval(['spiketimes_thiscell = C1.unitAC' num2str(u) ';']);
                
                ISI = log10( diff(spiketimes_thiscell) );
                ISI_piled{ct} = ISI;
                
                [N,edges] = histcounts(ISI, linspace(-4,3,200));
                Counts{ct} = N;
                Edges{ct} = N;
                
                E = edges(1:end-1);
                ind = find( E < - 2.5 );
                ind2 = find( E > -2.5);
                p{ct} = sum( N(ind) )./sum(N(ind2) )
                
                
                [fitobject, gof] = fit( edges(1:end-1)',N','gauss2');
                
                R2{ct} = gof.rsquare;
                NormW1{ct} = fitobject.a1./(fitobject.a1+fitobject.a2);
                NormW2{ct} = fitobject.a2./(fitobject.a1+fitobject.a2);
            end
        end;
    catch err
    end
end;





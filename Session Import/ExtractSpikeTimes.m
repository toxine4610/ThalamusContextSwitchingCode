
clear all;
block = 0;


%%
ValidNum = [];

X = dir('*.mat');
for i = 1:size(X,1)
    name = sprintf( X(i).name );
    target1 = strfind(name, 'T');
    target2 = strfind(name, '.');
    num = name( target1+2 : target2-1  );
    ValidNum = [ValidNum, str2num(num)];
end;

ValidNum = sort(ValidNum);

%%



for T = ValidNum
    
    warning off;
    D = '2018-07-27';
    mouseID = 'BPM43';
    format  = 'ExtHDD';
    
    if strcmpi(D, '2018-01-08') == 1; block = 0; end;
    
    disp(['Extracting data from --> TT ' num2str(T)] );
    
    if block == 0
        folder1 = 'X:\Rajeev\Cheetah\M4-3\2018-07-27_15-16-26\GlobalpeakC11\PFCMDSection1TT';
        C1 = load([folder1 num2str(T) '.mat']);
    elseif block == 1
        folder1 = '/Volumes/Data/somcre/2018-01-07_12-12-41/GlobalpeakC1/07-Jan-2018_PFCMDSection1TT';
        C1 = load([folder1 num2str(T) '.mat']);
        folder2 = '/Volumes/Data/somcre/2018-01-07_12-12-41/GlobalpeakC2/07-Jan-2018_PFCMDSection1TT';
        C2 = load([folder2 num2str(T) '.mat']);
    end;
    
    switch mouseID
        case 'SOMCRE1'
            PFCTT = 1:8;
            MDTT  = 9:16;
        case 'BPM33'
            PFCTT = [16:23];
            MDTT  = [1:15, 24];
    end;
    
    Z_C1 = C1.D;
   if block==1; Z_C2 = C2.D; end;
    
    num_units = size(C1.num_seq,1);
    
    if exist('Spfc')
        ct1 = numel(Spfc);
    else
        ct1 = 0;
    end;
    if exist('Smd')
        ct2 = numel(Smd);
    else
        ct2 = 0;
    end;
    
    for n = 1:num_units
        
        tetrode_num = C1.num_seq( n,1 );
        
        if ismember( tetrode_num ,PFCTT)
            ct1 = ct1+1;
            
            Spfc(ct1).filename = ['\' mouseID '\' D];
            Spfc(ct1).unitNbr = n;
            Spfc(ct1).ttNbr   = tetrode_num;
            
            if block ==0 || block==1
                % === process Par11 : correct
                
                eval( ['foo = C1.cl' num2str(n) '_cellDatabaseline;'] );
                clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpikeTimes_AC
                if isfield(foo, 'Vis_corr')
                    SpikeTimes_VC = foo.Vis_corr;
                else
                    SpikeTimes_VC = [];
                end;
                
                if isfield(foo, 'VisMix_corr')
                    SpikeTimes_VMC = foo.VisMix_corr;
                else
                    SpikeTimes_VMC = [];
                end;
                
                if isfield(foo, 'Aud_corr')
                    SpikeTimes_AC = foo.Aud_corr;
                else
                    SpikeTimes_AC = [];
                end;
                
                if isfield(foo, 'AudMix_corr')
                    SpikeTimes_AMC = foo.AudMix_corr;
                else
                    SpikeTimes_AMC = [];
                end;
                
                if block == 1
                    Spfc(ct1).SpikeTimes_Mix_Audition = [SpikeTimes_AC, SpikeTimes_AMC];
                    Spfc(ct1).SpikeTimes_Mix_Vision   = [SpikeTimes_VC, SpikeTimes_VMC];
                elseif block == 0
                    Spfc(ct1).SpikeTimes_Mix_Audition = [SpikeTimes_AMC];
                    Spfc(ct1).SpikeTimes_Mix_Vision   = [SpikeTimes_VMC];
                    Spfc(ct1).SpikeTimes_Pure_Audition = [SpikeTimes_AC];
                    Spfc(ct1).SpikeTimes_Pure_Vision   = [SpikeTimes_VC];
                end;
                
                clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
                
            end;
            
            if block == 1
                % === process Par12 : correct
                eval( ['foo = C2.cl' num2str(n) '_cellDatabaseline;'] );
                clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
                if isfield(foo, 'Vis_corr')
                    SpikeTimes_VC = foo.Vis_corr;
                else
                    SpikeTimes_VC = [];
                end;
                
                if isfield(foo, 'VisMix_corr')
                    SpikeTimes_VMC = foo.VisMix_corr;
                else
                    SpikeTimes_VMC = [];
                end;
                
                if isfield(foo, 'Aud_corr')
                    SpikeTimes_AC = foo.Aud_corr;
                else
                    SpikeTimes_AC = [];
                end;
                
                if isfield(foo, 'AudMix_corr')
                    SpikeTimes_AMC = foo.AudMix_corr;
                else
                    SpikeTimes_AMC = [];
                end;
                
                if block == 1
                    Spfc(ct1).SpikeTimes_Pure_Audition = [SpikeTimes_AC, SpikeTimes_AMC];
                    Spfc(ct1).SpikeTimes_Pure_Vision   = [SpikeTimes_VC, SpikeTimes_VMC];
                end;
                clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
            end;
            
            
        elseif ismember( tetrode_num, MDTT)
            ct2 = ct2+1;
            
            Smd(ct2).filename = ['\' mouseID '\' D];
            Smd(ct2).unitNbr = n;
            Smd(ct2).ttNbr   = tetrode_num;
            
            
            if block == 0 || block == 1
                % === process Part1 : correct
                
                eval( ['foo = C1.cl' num2str(n) '_cellDatabaseline;'] );
                clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
                if isfield(foo, 'Vis_corr')
                    SpikeTimes_VC = foo.Vis_corr;
                else
                    SpikeTimes_VC = [];
                end;
                
                if isfield(foo, 'VisMix_corr')
                    SpikeTimes_VMC = foo.VisMix_corr;
                else
                    SpikeTimes_VMC = [];
                end;
                
                if isfield(foo, 'Aud_corr')
                    SpikeTimes_AC = foo.Aud_corr;
                else
                    SpikeTimes_AC = [];
                end;
                
                if isfield(foo, 'AudMix_corr')
                    SpikeTimes_AMC = foo.AudMix_corr;
                else
                    SpikeTimes_AMC = [];
                end;
                
                if block == 1
                    Smd(ct2).SpikeTimes_Mix_Audition = [SpikeTimes_AC, SpikeTimes_AMC];
                    Smd(ct2).SpikeTimes_Mix_Vision   = [SpikeTimes_VC, SpikeTimes_VMC];
                elseif block == 0
                    Smd(ct2).SpikeTimes_Mix_Audition = [SpikeTimes_AMC];
                    Smd(ct2).SpikeTimes_Mix_Vision   = [SpikeTimes_VMC];
                    Smd(ct2).SpikeTimes_Pure_Audition = [SpikeTimes_AC];
                    Smd(ct2).SpikeTimes_Pure_Vision   = [SpikeTimes_VC];
                end;
                clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
                
            end;
            
            if block == 1
                % === process Part2 : correct
                eval( ['foo = C2.cl' num2str(n) '_cellDatabaseline;'] );
                clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
                if isfield(foo, 'Vis_corr')
                    SpikeTimes_VC = foo.Vis_corr;
                else
                    SpikeTimes_VC = [];
                end;
                
                if isfield(foo, 'VisMix_corr')
                    SpikeTimes_VMC = foo.VisMix_corr;
                else
                    SpikeTimes_VMC = [];
                end;
                
                if isfield(foo, 'Aud_corr')
                    SpikeTimes_AC = foo.Aud_corr;
                else
                    SpikeTimes_AC = [];
                end;
                
                if isfield(foo, 'AudMix_corr')
                    SpikeTimes_AMC = foo.AudMix_corr;
                else
                    SpikeTimes_AMC = [];
                end;
                
                if block == 1
                    Smd(ct2).SpikeTimes_Pure_Audition = [SpikeTimes_AC, SpikeTimes_AMC];
                    Smd(ct2).SpikeTimes_Pure_Vision   = [SpikeTimes_VC, SpikeTimes_VMC];
                end;
                clear SpikeTimes_AMC SpikeTimes_VMC SpikeTimes_VC SpiclkeTimes_AC
            end;
            
        end;
        
    end;
    
    clearvars -except Spfc Smd D Z_C1 Z_C2 block
    
end

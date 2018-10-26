
%%


SaveDir = [ cd '\' 'Randomized' '_Decode\'];

if exist( [SaveDir 'MD'] ) ~= 0
    rmdir([SaveDir 'MD'],'s')
    mkdir([SaveDir 'MD']);
elseif exist( [SaveDir 'MD'] ) == 0
    mkdir([SaveDir 'MD']);
end

if exist( [SaveDir 'PFC'] ) ~= 0
    rmdir([SaveDir 'PFC'],'s')
    mkdir([SaveDir 'PFC']);
elseif exist( [SaveDir 'PFC'] ) == 0
    mkdir([SaveDir 'PFC']);
end

%%

datapath = '../ForContextSwitchProject/DoubleCueDatabase/Production Randomized Trials';

X = dir(datapath);
MouseID = {X(3:end).name};
MyMouse  = { 'BPM13', 'BPM25', 'BPM33', 'BPM43' };

SmdMega = [];
SpfcMega = [];

ct = 0;

for m = 1:length(MouseID)
    mouse = MouseID{m};
    fprintf('Processing Mouse: %s\n', mouse );
    dates_this_mouse = getDateList([datapath filesep mouse]);
    for d = 1:length(dates_this_mouse)
        fprintf('  Processing date: %s\n', dates_this_mouse{d});
        ST = load( [datapath filesep mouse filesep dates_this_mouse{d} filesep MyMouse{m} '_randomized_session' ] );
        
        clear Spfc_mix Smd_mix ZVis ZAud
        
        [Spfc_mix, Smd_mix, ZVis, ZAud] =  packageData(ST.Z_C1, ST.Smd, ST.Spfc);
        
        % ==== separate into trial type.
        for i = 1:numel(Smd_mix)
            
            clear raster_data raster_labels raster_site_info
            keep = 1;
            if keep == 1
                [raster_data,raster_labels,raster_site_info] = Mixed_prepDateforDecodingToolbox(Smd_mix(i),d);
                save( [SaveDir 'MD\MD_sess_' num2str(d) '_cell_' num2str(i) '_Raster_Data'], 'raster_data','raster_labels','raster_site_info');
            elseif keep == 0
                rejectMD = [rejectMD, i];
            end;
        end;
        
        for i = 1:numel(Spfc_mix)
            clear raster_data raster_labels raster_site_info
            keep = 1;
            if keep == 1
                [raster_data,raster_labels,raster_site_info] = Mixed_prepDateforDecodingToolbox(Spfc_mix(i),d);
                save( [SaveDir 'PFC\PFC_sess_' num2str(d) '_cell_' num2str(i) '_Raster_Data'], 'raster_data','raster_labels','raster_site_info');
            end;
        end
    end
end

%% BIN Data........................................................

% add the path to the Neural Decoding Toolbox
toolbox_directory_name = 'C:\Users\Halassalab-CG\Dropbox (Personal)\Rajeev\Rajeev_Code\ndt.1.0.4_exported\';
addpath(toolbox_directory_name);
add_ndt_paths_and_init_rand_generator;


bin_width = 100;
step_size = 10;
decodeType = 'context';


switch decodeType
    case 'Context'
        dc = 'trial_context';
    case 'TrialType'
        dc = 'stimulus_ID';
    case 'RuleMeaning'
        dc = 'rule_meaning';
    case 'all'
        dc = 'cominedInfo';
end;

%%

resultsFolder = [ cd '\results_MD_DC_' date '_\'];
mkdir(resultsFolder);
mkdir( [resultsFolder '\shuff_results'] );

raster_data_directory_name  = [SaveDir, 'MD\'];
save_prefix_name = ['Binned_MD_all'];

binned_data_file_name = create_binned_data_from_raster_data(raster_data_directory_name, save_prefix_name, bin_width, step_size);


numUnits =  length( dir( [raster_data_directory_name '*.mat'] ) );
units_to_sample = numUnits;


%%
for s = 1:length(units_to_sample)
    fprintf('Random Sampling Size = %d (of % d MD neurons)\n', units_to_sample(s), numUnits);
    for iter = 1:20
        u = randsample(1:numUnits, units_to_sample(s));
        DECODING_RESULTS  = Context_run_basic_decoding_shuff(0, binned_data_file_name, 'MD', dc,  u, resultsFolder);
        MDMI{s}(iter,:) = DECODING_RESULTS.MUTUAL_INFORMATION.from_combined_confusion_matrix_over_all_resamples.decoding_results;
        MDCA{s}(iter,:) = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results;
        clear DECODING_RESULTS
    end;
end;

%% PFC ===================================================================

resultsFolder = [ cd '\results_PFC_DC_' date '_\'];
mkdir(resultsFolder);
mkdir( [resultsFolder '\shuff_results'] );

raster_data_directory_name  = [SaveDir, 'PFC\'];
save_prefix_name = ['Binned_PFC_all'];

binned_data_file_name = create_binned_data_from_raster_data(raster_data_directory_name, save_prefix_name, bin_width, step_size);

numUnits =  length( dir( [raster_data_directory_name '*.mat'] ) );
units_to_sample = numUnits;


for s = 1:length(units_to_sample)
    fprintf('Random Sampling Size = %d (of % d PFC neurons)\n', units_to_sample(s), numUnits);
    for iter = 1:20
        u = randsample(1:numUnits, units_to_sample(s));
        DECODING_RESULTS  = Context_run_basic_decoding_shuff(0, binned_data_file_name, 'PFC', dc,  u, resultsFolder);
        PFCMI{s}(iter,:) = DECODING_RESULTS.MUTUAL_INFORMATION.from_combined_confusion_matrix_over_all_resamples.decoding_results;
        PFCCA{s}(iter,:) = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results;
        clear DECODING_RESULTS
    end;
end;


%%
colorMD = [0.49,0.18,0.56];
colorPFC = [0.5,0.5,0.5];

t = linspace(-0.2,1.2, length(MDMI{1}(1,:)) ) ;
clrmap = lines(5);
figure(2); set(gcf,'color','w')

for i = 1:length(units_to_sample)
    stdshadenew( t, mean( PFCCA{i} ), std( PFCCA{i} ), 0.5, colorPFC);
    stdshadenew( t, mean( MDCA{i} ), std( MDCA{i} ), 0.5, colorMD);
    axis square; box off; set(gca,'tickdir','out','fontsize',16);
end;
xlim([-0.2, 1.2] );
set(gcf, 'position', [250 300 950 300])  % expand the figure;


%%

ind1= find(t >= -0.2 & t <= 0.9);
ind2= find(t >= -0.2 & t <= 1);
for i = 1:20
    aucMD(i) = trapz( MDCA{1}(i,ind1) );
    aucPFC(i) = trapz( PFCCA{1}(i,ind2));
end;

figure(3); set(gcf,'color','w');
bh = boxplot([aucMD; aucPFC]' );

set( bh(1:4,2),'color', colorPFC,'linewidth',3);
set( bh(1:4,1),'color', colorMD,'linewidth',3);
set( bh(6,:),'color', 'k','linewidth',3);

axis square; box off; set(gca,'tickdir','out','fontsize',16)

set( bh(5,1),'color', colorMD,'linewidth',3);
set( bh(5,2),'color', colorPFC,'linewidth',3);

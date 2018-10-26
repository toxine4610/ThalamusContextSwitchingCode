the_classifier = max_correlation_coefficient_CL;


% the_feature_preprocessors{1} = zscore_normalize_FP;


resultsFolder = [ cd '\results_MD_DC_' date '_\'];
mkdir(resultsFolder);
mkdir( [resultsFolder '\shuff_results'] );

raster_data_directory_name  = [SaveDir, 'MD'];
save_prefix_name = ['Binned_MD_all'];

binned_data_file_name = create_binned_data_from_raster_data(raster_data_directory_name, save_prefix_name, bin_width, step_size);


numUnits =  length( dir( [raster_data_directory_name '*.mat'] ) );
units_to_sample = 200;

num_cv_splits =5;

id_string_names = {'BlueNoise', 'BrownNoise', 'BlueLED', 'GreenLED'};
% pos_string_names = {'BlueNoise', 'Brown_Noise'};


for iTrainPosition = 1:4
    for iTestPosition = 1:4
        
        
        for iID = 1:4
%             the_training_label_names{iID} = {[id_string_names{iID} '_' pos_string_names{iTrainPosition}]};
%              the_test_label_names{iID} =  {[id_string_names{iID} '_' pos_string_names{iTestPosition}]};
            
            the_training_label_names{iID} = {[id_string_names{iID}]};
            the_test_label_names{iID} =  {[id_string_names{iID}]};
        end
        
        ds = generalization_DS(binned_data_file_name, 'stimulus_ID', num_cv_splits, the_training_label_names, the_test_label_names);
        the_cross_validator.test_only_at_training_times = 1;

        the_cross_validator = standard_resample_CV(ds, the_classifier);
        the_cross_validator.num_resample_runs =2;
        DECODING_RESULTS = the_cross_validator.run_cv_decoding;
        all_results{iTrainPosition, iTestPosition} = DECODING_RESULTS;
        
    end
end

%%
t = linspace(-0.2,1.2,131);
ind = find( t>=0 & t<0.8);
indB = find( t>=-1 & t < -0.5);

figure(1); set(gcf,'color','w');

position_names = {'BLNoise', 'BRNoise','BlueLED','GreenLED'};

for iTrainPosition = 1:4
    % load the results from each training and test location
    for iTestPosition = 1:4
        AR(iTrainPosition, iTestPosition) = max( all_results{iTrainPosition,iTestPosition}.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(ind) );
    end
    
    % create a bar plot for each training lcoation
    subplot(1, 4,  iTrainPosition)
    bar(AR(iTrainPosition, :) .* 100);
    
    title(['Train ' position_names{iTrainPosition}])
    ylabel('Classification Accuracy');
    set(gca, 'XTickLabel', position_names)
    xlabel('Test position')
    xLims = get(gca, 'XLim');
    line([xLims], [1/2 1/2] .* 100, 'color', [0 0 0]);    % put a line at the chance level of decoding
end;
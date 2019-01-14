
% This script executes plotting code for figure 5 panels g,h of Yates et al. (2017)
% To run successfully, you must run cmdFitLIPtrunc.m first



% path to the data base directory. dataPath should have subdirectories
% neurons, stim, main_fits, lip_trunc_fits
dataPath = getpref('mtlipglm', 'dataPath');

%% Load up LIP fits

% The raw data have been fit with the required models and analyzed. To
% regenerate fits, use cmdFitLIP.m
% Note: The fits in cmdFitLIP.m are generated with a fixed hyperparameter 
% so the exact numbers may subtly differ from the published text

% P = loadUpFits(dataPath, 'LIP', 'lip_trunc_fits');
P = loadUpFits(dataPath, 'LIP', 'glmLIPtruncatedChoice');



%% Plot kernels and goodnes-of-fit
figure(1); clf
figure(2); clf


nModels=numel(P(1).model)-1;

% build color map
aa=10;
cmap1=(cbrewer('seq', 'Reds', nModels+aa));
cmap2=(cbrewer('seq', 'Blues', nModels+aa));

cmap1(1:aa,:)=[];
cmap2(1:aa,:)=[];

r2    =nan(nModels,1);
r2SD  =nan(nModels,1);
truncations=nan(nModels,1);

% exponentiate kernels, or not
f=@exp; %f=@(x) x;
for kModel=1:nModels
    
    % plot kernels
    figure(1)
    
    field = 'resp1';
    kernels = cell2mat(arrayfun(@(x) mean(x.model(kModel+1).wts.(field).data,2)'*-sign(x.model(1).dprime), P, 'UniformOutput', false)');
    tr = P(1).model(kModel+1).wts.(field).tr;
    plot(P(1).model(kModel+1).wts.(field).tr, mean(f(kernels)), 'Color', cmap1(kModel,:)); hold on
    
    field='resp2';
    kernels=cell2mat(arrayfun(@(x) mean(x.model(kModel+1).wts.(field).data,2)'*-sign(x.model(1).dprime), P, 'UniformOutput', false)');
    plot(tr, mean(f(kernels)), 'Color', cmap2(kModel,:)); hold on
    
    % plot goodness-of-fit
    figure(2)    
    r2(kModel)   = mean(arrayfun(@(x) x.model(kModel+1).r2psth, P));
    r2SD(kModel) = std(arrayfun(@(x) x.model(kModel+1).r2psth, P)) / sqrt(numel(P));
    
    truncations(kModel)=min(tr);
    
    errorbar(min(tr), r2(kModel), r2SD(kModel), '.', 'Color', cmap2(kModel,:)); hold on
end




figure(1); 
xlabel('Time (aligned to saccade)')
ylabel('weight')
xlim([-2500 500])

figure(2);
plot(truncations, r2, 'k')
xlabel('Time included before saccade')
ylabel('% variance explained PSTH')
xlim([-2500 0])

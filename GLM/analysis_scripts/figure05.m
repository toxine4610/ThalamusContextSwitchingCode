% This script executes plotting code for figure 5 panels b-f of Yates et al. (2017)
% To plot panels g,h, use figure05trunc.m

% path to the data base directory. dataPath should have subdirectories
% neurons, stim, main_fits, lip_trunc_fits
dataPath = getpref('mtlipglm', 'dataPath');

%% Load up LIP fits

% The raw data have been fit with the required models and analyzed. To
% regenerate fits, use cmdFitLIP.m
% Note: The fits in cmdFitLIP.m are generated with a fixed hyperparameter 
% so the exact numbers may subtly differ from the published text

S = loadUpFits(dataPath, 'LIP');

%% Plot Model PSTH on top of data PSTH
% The three models compared in figure 5 are the "Stimulus-to-LIP" model,
% the "MT-to-LIP" model and "MT-to-LIP (full choice)"

modelNames=arrayfun(@(x) x.name, S(1).model, 'UniformOutput', false);

% The three models have different names in the code:
% Stimulus-to-LIP = Poisson
% MT-to-LIP       = MTsim
% MT-to-LIP (full choice) = MTsimChoice
% Compare them to the data in each plot
modelIxs={[find(strcmp(modelNames, 'data')) find(strcmp(modelNames, 'Poisson'))], ...
    [find(strcmp(modelNames, 'data')) find(strcmp(modelNames, 'MTsim'))], ...
    [find(strcmp(modelNames, 'data')) find(strcmp(modelNames, 'MTsimChoice'))]};


for iModel=1:numel(modelIxs)
    figure(iModel); clf
    
    modelIx  = modelIxs{iModel};
    nModels  = numel(modelIx);
    nNeurons = numel(S);
    
    plotFields={'psthCoh', 'ptaRaw'};
    nPlotFields = numel(plotFields);
    
    for f=1:nPlotFields
        subplot(1,nPlotFields, f)
        
        field=plotFields{f};
        
        
        % --- Change Color Map
        if any(strfind(field, 'pta'))
            cmap=hot(12);
        else
            cmap=cbrewer('jake', 'rdbu', 8);
        end
        
        for kM=1:nModels
            
            kModel=modelIx(kM);
            
            % get normalized population response
            cohpData = reshape(cell2mat(arrayfun(@(x) x.model(kModel).(field)/max(x.model(1).(field)(:)), S, 'UniformOutput', false)), [size(S(1).model(kModel).(field)) nNeurons]);
            popcoh   = nanmean(cohpData,3);
            
            % set time axis
            if any(strfind(field, 'pta'))
                timex=S(1).model(1).ptaTime;
                xdim = [0 2000];
            else
                timex=S(1).model(1).psthTime(:);
                xdim = [-500 1500];
            end
            
            % change marker for data and model
            if kModel==1
                marker='.';
            else
                marker='-';
            end
            
            % plot
            for kP=1:size(popcoh,2)
                if size(timex,2)>1
                    txx=timex(:,kP);
                else
                    txx=timex;
                end
                
                plot(txx, popcoh(:,kP), marker, 'Color', cmap(kP,:), 'MarkerSize', 2); hold on
            end
            
            % title is model name
            title(S(1).model(kModel).name)
            xlim(xdim)
            
        end
        
    end
    
end

%% Compare models Goodness-of-fit

% --- variance explained of the PSTH
field='r2psth';
modelNames={S(1).model.name};
r2MT    = arrayfun(@(x) x.model(strcmp(modelNames, 'MTsim')).(field), S);
r2St    = arrayfun(@(x) x.model(strcmp(modelNames, 'Poisson')).(field), S);
r2MTcho = arrayfun(@(x) x.model(strcmp(modelNames, 'MTsimChoice')).(field), S);

figure(10); clf
subplot(1,2,1)
plot(r2St, r2MT, '.', xlim, xlim, '-')
xlabel('Stimulus-To-LIP')
ylabel('MT-To-LIP')

subplot(1,2,2)
plot(r2MT, r2MTcho, '.', xlim, xlim, '-')
xlabel('MT-To-LIP')
ylabel('MT-To-LIP (Full Choice Term)')

% --- test log likelihood 
field='llpoisson';
llMT    = arrayfun(@(x) x.model(strcmp(modelNames, 'MTsim')).(field), S);
llSt    = arrayfun(@(x) x.model(strcmp(modelNames, 'Poisson')).(field), S);
llMTcho = arrayfun(@(x) x.model(strcmp(modelNames, 'MTsimChoice')).(field), S);

figure(11); clf
subplot(1,2,1)
histogram(llMT-llSt);
xlabel('MT-To-LIP / Stimulus-To-LIP')

subplot(1,2,2)
histogram(llMTcho-llMT);
xlabel('Full Choice Term / Truncated')

[stest, ~, stats]=signrank(r2MT-r2St);
fprintf(1, 'MT input better than Stim by %02.0f%%, sign test (zval=%d, sign=%d) p=%d\n', nanmean(r2MT-r2St)*100, stats.zval, stats.signedrank, stest);
[stest, ~, stats]=signrank(r2MTcho-r2MT);
fprintf(1, 'Choice better than no Choice by %02.0f%%, sign test (zval=%d, sign=%d) p=%d\n', nanmean(r2MTcho-r2MT)*100, stats.zval, stats.signedrank, stest);
fprintf(1, '%02.2f %% +- %02.2f of variance of PSTH explained by MTcho\n', 100*mean(r2MTcho), 100*std(r2MTcho)/sqrt(nNeurons));
[stest, ~, stats]=signrank(llMT-llSt);
fprintf(1, 'MT (log-likelihood) input better than Stim by %02.3f, sign test (zval=%d, sign=%d) p=%d\n', mean(llMT-llSt), stats.zval, stats.signedrank, stest);


%% Plot Kernels
plotAll=false;

figure(12); clf

% --- Targets
subplot(3,1,2)
kModel=strcmp(modelNames, 'MTsimChoice');
field='targson';
f=@exp; % exponentiate for gain
kernels=cell2mat(arrayfun(@(x) mean(x.model(kModel).wts.(field).data,2)'*sign(sum(mean(x.model(kModel).wts.(field).data,2))), S, 'UniformOutput', false)');

m=mean(f(kernels));

tx=S(1).model(kModel).wts.(field).tr;
if plotAll
    plot(tx, f(kernels), 'Color', .5*[1 1 1]); hold on
end
plot(tx, m, 'k', 'Linewidth', 1); hold on
xd=[0 1500];
plot(xd, f([0 0]), 'k', 'Linewidth', .5)

title('Targets')
ylabel('Gain')

% --- Choice
subplot(3,1,3)
field='resp1';
Cho1=cell2mat(arrayfun(@(x) mean(x.model(kModel).wts.(field).data,2)', S, 'UniformOutput', false)')';
field='resp2';
Cho2=cell2mat(arrayfun(@(x) mean(x.model(kModel).wts.(field).data,2)', S, 'UniformOutput', false)')';

% use dprime measure to align to neuron's preferred direction
dp = arrayfun(@(x) x.model(1).dprime, S);
prefDir=-sign(dp); %sum(Cho2-Cho1)); % computed from preferred choice direction
InRF=bsxfun(@times, Cho2, prefDir); % align to the preferred direction
OutRF=bsxfun(@times, Cho1, prefDir);

InRF=f(InRF);
OutRF=f(OutRF);

m=mean(InRF,2);

tx=S(1).model(kModel).wts.(field).tr;
cmap=flipud(cbrewer('jake', 'rdbu', 2));
if plotAll
    for i=1:size(InRF,2)
        plot(tx, InRF(:,i),  'Linewidth', .1, 'Color', cmap(1,:)); hold on
        plot(tx, OutRF(:,i), 'Linewidth', .1, 'Color', cmap(2,:));
    end
else
plot(tx, m, 'Color', cmap(1,:), 'Linewidth', 1); hold on
m=mean(OutRF,2);
plot(tx, m, 'Color', cmap(2,:), 'Linewidth', 1);
end

%------------
% plot truncated choice kernels
subplot(3,1,3)
kModel=strcmp(modelNames, 'MTsim');

field='resp1';
Cho1=cell2mat(arrayfun(@(x) mean(x.model(kModel).wts.(field).data,2)', S, 'UniformOutput', false)')';
field='resp2';
Cho2=cell2mat(arrayfun(@(x) mean(x.model(kModel).wts.(field).data,2)', S, 'UniformOutput', false)')';

InRF=bsxfun(@times, Cho2, prefDir); % align to the preferred direction
OutRF=bsxfun(@times, Cho1, prefDir);

InRF=f(InRF);
OutRF=f(OutRF);

if ~plotAll
m=mean(InRF,2);
tx=S(1).model(kModel).wts.(field).tr;
plot(tx, m, ':', 'Color', cmap(1,:), 'Linewidth', .5); hold on
m=mean(OutRF,2);
plot(tx, m, ':', 'Color', cmap(2,:), 'Linewidth', .5);
end
% save with axes
axis tight
if f(.5)==1
    
if plotAll
    ylim([-2.5 3.5])
else
    ylim([-.4 .8])
end
set(gca, 'YTick', min(ylim):.5:max(ylim))
else
   ylim([.5 3]) 
end
xlim([-2500 500])
xlabel('Time from covariate onset (ms)')


title('Choice')

% -------------------------------------------------------------------------
% --- MT kernels
% This code assesses the magnitude of a single pulse of simulated
% MT input on LIP spike rate in units of spike rate gain

subplot(3,1,1)

% --- First, load the stimulus kernels from the MT fits
load(fullfile(dataPath, 'MTkernels.mat'))

% kCon = mean(kCons(:,isPat),2);
% kDir = mean(kDirs(:,isPat),2);
% b0 = mean(b0(isPat));

% take the average as an approximation to an idealized MT neuron
kCon = mean(kCons,2);
kDir = mean(kDirs,2);
b0 = mean(b0);

% check how the fits were generated and adjust stimulus accordingly
if contrastIsBoxcar
    contrast=[zeros(1,20), ones(1,105), zeros(1,20)];
else
    contrast=[zeros(1,20), .2*ones(1,105), zeros(1,20)];
end

% simulated direction input
pulse=[zeros(1, 20) (1/19)*ones(1,15) zeros(1, 20 + 6*15)];

% MT is simulated by convolving the filters with the stimulus input,
% summing and then exponentiating
Xdir1 = filter(kDir, 1, -pulse);
Xdir0 = filter(kDir, 1, pulse);
Xcon = filter(kCon, 1, contrast);
% pulse in the preferred direction

binSize = diff(S(1).model(2).psthTime(1:2));
MTprefPulse1 = exp(b0 + Xdir1 + Xcon)/binSize;
MTantiPulse1 = exp(b0 -Xdir1 + Xcon)/binSize;
% pulse in the antipreferred direction
MTprefPulse0 = exp(b0 + Xdir0 + Xcon)/binSize;
MTantiPulse0 = exp(b0 -Xdir0 + Xcon)/binSize;

% Load MT-to-LIP filters from model fits
kModel=strcmp(modelNames, 'MTsim');
field='MTpref';
MTpref=cell2mat(arrayfun(@(x) mean(x.model(kModel).wts.(field).data,2)', S, 'UniformOutput', false)')';
field='MTanti';
MTanti=cell2mat(arrayfun(@(x) mean(x.model(kModel).wts.(field).data,2)', S, 'UniformOutput', false)')';

% correct for the preferred choice of each LIP unit
MTpref=bsxfun(@times, MTpref, prefDir); % align to the preferred direction
MTanti=bsxfun(@times, MTanti, prefDir);

% convolve MT-to-LIP filters with simulated MT input
Cpref = convmtx(MTprefPulse1(:), size(MTpref, 1));
Canti = convmtx(MTantiPulse1(:), size(MTpref, 1));

% convolve with MT input
MTpref1 = Cpref*MTpref;
MTanti1 = Canti*MTanti;

Cpref = convmtx(MTprefPulse0(:), size(MTpref, 1));
Canti = convmtx(MTantiPulse0(:), size(MTpref, 1));
MTpref0 = Cpref*MTpref;
MTanti0 = Canti*MTanti;

tx = (1:size(MTpref1,1))*binSize - 200;
% Effect of a pulse in the preffered direction over a pulse in the
% anti-preffered direction
MTpulseEffect=f(MTpref1+MTanti1 - MTpref0-MTanti0);
plot(tx, mean(MTpulseEffect,2)); hold on

ylim([.95 1.1])
xlim([0 1.5e3])
set(gca, 'YTick', min(ylim):.05:max(ylim), 'XTick', 0:250:1.5e3)

title('simulated MT input')




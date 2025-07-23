function DecodeMulticlass(cfg0,subject)

%% Save Path
outputDir = fullfile(cfg0.root,cfg0.output_path,subject);
if ~exist(outputDir,'dir'); mkdir(outputDir); end


%% Load MEG data
disp('loading..')
disp(subject)
if contains(cfg0.task,'sym')
    data = load(fullfile(cfg0.root,'EpochedData',subject,'sym_trials.mat'));
    data = data.sym_trials;
    time = data.time{1};
    disp('loaded data')
elseif contains(cfg0.task,'dot')
    data = load(fullfile(cfg0.root,'EpochedData',subject,'dot_trials.mat'));
    data = data.dot_trials;
    time = data.time{1};
    disp('loaded data')
elseif contains(cfg0.task,'det')
    data = load(fullfile(cfg0.root,'EpochedData',subject,'det_trials.mat'));
    data = data.det_trials;
    time = data.time{1};
    disp('loaded data')
end

%% Remove No Resp Trials
cfgS = [];
cfgS.trials = data.trialinfo(:,3) ~= 0 ;
data = ft_selectdata(cfgS,data);

%% Get Trial x Channels x Time Matrix
cfgS = [];
cfgS.keeptrials = true;
cfgS.channel=cfg0.channel;
data = ft_timelockanalysis(cfgS,data);

%% Select data
if isfield(cfg0,'selData')
    cfgS = [];
    cfgS.trials = eval(cfg0.selData);
    data = ft_selectdata(cfgS,data);
end

%% Smooth Data
smoothed_data = zeros(size(data.trial));
for trial = 1:size(data.trial,1)
    smoothed_data(trial,:,:) = ft_preproc_smooth(squeeze(data.trial(trial,:,:)),cfg0.nMeanS);
end

%% Get labels
labels = eval(cfg0.label);
labels = labels - min(labels) + 1; %make sure labels start at 1

cfgS = [];
cfgS.classifier = 'multiclass_lda';
cfgS.metric = cfg0.metric;
cfgS.preprocess ={'undersample','average_samples'};
cfgS.preprocess_param{2}.group_size = cfg0.nAvgSamples;
cfgS.repeat = 1;
cfgS.stratify = true;
[results,~] = mv_classify_timextime(cfgS,smoothed_data,labels);

%% Save
save(fullfile(outputDir,[cfg0.output_prefix '_' cfg0.metric]),'results');


%% Plot
if cfg0.plot
    % Calculate color limits
    center_value = 1/6;  % Desired center point for color transition
    c_min = min([min(results,[],'all'),min(results,[],'all')])-0.05;
    c_max = max([max(results,[],'all'),max(results,[],'all')])+0.05;

    % Make the colormap symmetric around the center value
    c_range = max(abs([c_min-center_value, c_max-center_value]));
    caxis_range = [center_value-c_range, center_value+c_range];

    figure('units','normalized','outerposition',[0 0 1 1])

    subplot(2,1,1)
    wdw_size = 3;
    smoothresults = movmean(diag(results),wdw_size);
    plot(time, smoothresults, 'LineWidth', 1.5, 'Color', 'k')
    yline(1/unique(labels),'r--')
    xline(0,'r')
    xlabel('Time (s)'); ylabel('Accuracy');

    subplot(2,1,2)
    imagesc(time, time, results, caxis_range); axis xy; colorbar;
    xlabel('Train Time (s)'); ylabel('Test Time (s)');
    hold on; plot([0 0], ylim, 'r'); hold on; plot(xlim, [0 0], 'r');
    hold on; plot([1 1], ylim, 'k'); hold on; plot(ylim, [1 1], 'k');

    % Use flipped colormap and set color axis
    cmap = flipud(ft_colormap('RdBu'));
    colormap(cmap);
    caxis(caxis_range);  % Set the color limits with chance as the midpoint
    axis square
    sgtitle([subject,' ',cfg0.pltTitle])

    drawnow;

end
end
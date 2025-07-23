function CrossDecodeBinary(cfg0,subject)
%% Cross decode grating ocntrast with empty sets

%% Save Path
outputDir = fullfile(cfg0.root,cfg0.output_path,subject);
if ~exist(outputDir,'dir'); mkdir(outputDir); end 

%% Load MEG data
disp('loading..')
disp(subject)

dot_data = load(fullfile(cfg0.root,'EpochedData',subject,'dot_trials.mat'));
dot_data = dot_data.dot_trials;
time = dot_data.time{1};
disp('loaded dot data')

det_data = load(fullfile(cfg0.root,'EpochedData',subject,'det_trials.mat'));
det_data = det_data.det_trials;
time = det_data.time{1};
disp('loaded det data')


%% Remove No Resp Trials
cfgS = [];
cfgS.trials = det_data.trialinfo(:,3) ~= 0 ;
det_data = ft_selectdata(cfgS,det_data);
cfgS.trials = dot_data.trialinfo(:,3) ~= 0 ;
dot_data = ft_selectdata(cfgS,dot_data);


%% Get Trial x Channels x Time Matrix 
cfgS = [];
cfgS.keeptrials = true;
cfgS.channel=cfg0.channel;
det_data = ft_timelockanalysis(cfgS,det_data);
dot_data = ft_timelockanalysis(cfgS,dot_data);


%% Select data
if isfield(cfg0,'selData')
    cfgS = [];
    cfgS.trials = eval(cfg0.selData);
    det_data = ft_selectdata(cfgS,det_data);
end

%% Smooth Data
smoothed_dot_data = zeros(size(dot_data.trial));
smoothed_det_data = zeros(size(det_data.trial));

for trial = 1:size(dot_data.trial,1)
    smoothed_dot_data(trial,:,:) = ft_preproc_smooth(squeeze(dot_data.trial(trial,:,:)),cfg0.nMeanS);
end
for trial = 1:size(det_data.trial,1)
    smoothed_det_data(trial,:,:) = ft_preproc_smooth(squeeze(det_data.trial(trial,:,:)),cfg0.nMeanS);
end

%% Binarise labels
det_labels = eval(cfg0.det_label);
det_labels = det_labels - min(det_labels) + 1; %make sure labels start at 1
dot_labels = eval(cfg0.num_label);
dot_labels = dot_labels - min(dot_labels) + 1; %make sure labels start at 1

%% Undersample Non-Zero Class
%So each non-zero numerosity appears the same number of times
cfgB = [];
cfgB.numNTClass = 5;
[smoothed_dot_data,dot_labels] = UndersampleBinarise(cfgB,smoothed_dot_data, dot_labels, 1);

%% Turn all numbers to zero and not-zero
dot_labels = (dot_labels == 1)+1; %zero = 2
%det_labels = (det_labels == 2)+1; %absence = 2

%% Within Time x Time Decoding
cfgS = [];
cfgS.classifier = 'lda';
cfgS.metric = cfg0.metric;
cfgS.preprocess ={'undersample','average_samples'};
cfgS.preprocess_param{2}.group_size = cfg0.nAvgSamples;
cfgS.repeat = 1;
if isfield(cfg0,'balanceContrasts') & cfg0.balanceContrasts
    cfgS.cv = 'leaveout';
else
    cfgS.cv = 'kfold';
    cfgS.k  = 5;
    cfgS.stratify = true;
end
[results,~] = mv_classify_timextime(cfgS,smoothed_dot_data,dot_labels,smoothed_det_data,det_labels);
save(fullfile(outputDir,[cfg0.output_prefix '_' cfg0.metric '_trainDots']),'results');
[results,~] = mv_classify_timextime(cfgS,smoothed_det_data,det_labels,smoothed_dot_data,dot_labels);
save(fullfile(outputDir,[cfg0.output_prefix '_' cfg0.metric '_trainDet']),'results');


%% Plot
if cfg0.plot
    c_min = min([min(results,[],'all'),min(results,[],'all')])-0.05;
    c_max = max([max(results,[],'all'),max(results,[],'all')])+0.05;
    
    figure('units','normalized','outerposition',[0 0 1 1])

    subplot(2,1,1)
    wdw_size = 3;
    smoothresults = movmean(diag(results),wdw_size);
    plot(time,smoothresults,'LineWidth',1.5,'Color','k')
    yline(0.5,'r--')
    xline(0,'r')
    xlabel('Time (s)')
    ylabel('Accuracy')

    subplot(2,1,2)
    imagesc(time,time,results,[c_min,c_max]); axis xy; colorbar;
    xlabel('Train Time (s)'); ylabel('Test Time (s)');
    hold on; plot([0 0],ylim,'r'); hold on; plot(xlim,[0 0],'r');
    hold on; plot([1 1],ylim,'k'); hold on; plot(ylim,[1 1],'k')
    cmap = flipud(ft_colormap('RdBu'));
    colormap(cmap);
    axis square
    sgtitle([subject, ' ',cfg0.pltTitle])
    drawnow;
end



end
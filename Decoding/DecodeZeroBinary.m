function DecodeZeroBinary(cfg0,subject)

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

%% Binarise labels
labels = eval(cfg0.label);
labels = labels + 1; %make sure labels start at 1

%% Undersample Non-Zero Class
%So each non-zero numerosity appears the same number of times
cfgB = [];
cfgB.numNTClass = 5;
[smoothed_data,labels] = UndersampleBinarise(cfgB,smoothed_data, labels, 1);

%% Turn all numbers to zero and not-zero
labels = (labels == 1)+1; %zero = 2



%% Within Time x Time Decoding
cfgS = [];
cfgS.classifier = 'lda';
cfgS.metric = cfg0.metric;
cfgS.preprocess ={'undersample','average_samples'};
cfgS.preprocess_param{2}.group_size = cfg0.nAvgSamples;
cfgS.repeat = 1;
cfgS.cv = 'kfold';
cfgS.k  = 5;
cfgS.stratify = true;
[results,~] = mv_classify_timextime(cfgS,smoothed_data,labels);

%% Save
save(fullfile(outputDir,[cfg0.output_prefix '_' cfg0.metric]),'results');


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
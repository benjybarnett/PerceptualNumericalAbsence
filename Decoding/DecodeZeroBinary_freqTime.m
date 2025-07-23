function DecodeZeroBinary_freqTime(cfg0,subject,band)

%% Save Path
outputDir = fullfile(cfg0.root,cfg0.output_path,band,subject);
if ~exist(outputDir,'dir'); mkdir(outputDir); end 

%% Load MEG data
disp('loading..')
disp(subject)
if contains(cfg0.task,'sym')
    data = load(fullfile(cfg0.root,'EpochedBandData',band,subject,'sym_trials.mat'));
    data = data.sym_trials;
    time = data.time{1};
    load(fullfile(cfg0.root,'rejectedTrials',subject,'rej_trials_sym.mat'));
    disp('loaded data')
elseif contains(cfg0.task,'dot')
    data = load(fullfile(cfg0.root,'EpochedBandData',band,subject,'dot_trials.mat'));
    data = data.dot_trials;
    time = data.time{1};
    load(fullfile(cfg0.root,'rejectedTrials',subject,'rej_trials_dot.mat'));
    disp('loaded data')
elseif contains(cfg0.task,'det')
    data = load(fullfile(cfg0.root,'EpochedBandData',band,subject,'det_trials.mat'));
    data = data.det_trials;
    time = data.time{1};
    load(fullfile(cfg0.root,'rejectedTrials',subject,'rej_trials_det.mat'));

    disp('loaded data')
end

%% Remove the visually rejected artefacts
cfgS = [];
cfgS.trials = ~ismember(data.trialnumbers,rej_trials);
data = ft_selectdata(cfgS,data);


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


%% Hilbert Transform
[nTrials, nChans, nTime] = size(data.trial);

% Preallocate an array for the analytic signal:
X_alpha = zeros(size(data.trial));
ds = {data};
bands = {X_alpha};
% Loop over trials and channels
for b = 1:length(ds)
    band_data = ds{b}.trial;
    for t = 1:nTrials
        for c = 1:nChans
            % Extract the time series for the c-th channel in the t-th trial
            signal = squeeze(band_data(t, c, :));  % size will be (nTime,)
    
            % Take the Hilbert transform along time
            analyticSig = abs(hilbert(signal));
    
            % Store the analytic signal 
            bands{b}(t, c, :) = analyticSig;
        end
    end
end

%% Binarise labels
labels = eval(cfg0.label);
labels = labels + 1; %make sure labels start at 1

%% Undersample Non-Zero Class
%So each non-zero numerosity appears the same number of times
new_labels = {};
cfgB = [];
cfgB.numNTClass = 5;
for b = 1:length(bands)
    [bands{b},new_labels{b}] = UndersampleBinarise(cfgB,squeeze(bands{b}), labels, 1);
end


%% Turn all numbers to zero and not-zero
for i = 1:length(new_labels)
    new_labels{i} = (new_labels{i} == 1)+1; %zero = 2
end

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
[results{1},~] = mv_classify_timextime(cfgS,bands{1},new_labels{i});

%% Save
save(fullfile(outputDir,[cfg0.output_prefix '_' cfg0.metric]),'results');


%% Plot
if cfg0.plot
    c_min = min([min(results{1},[],'all'),min(results{1},[],'all')])-0.05;
    c_max = max([max(results{1},[],'all'),max(results{1},[],'all')])+0.05;
    
    figure('units','normalized','outerposition',[0 0 1 1])

    subplot(2,1,1)
    wdw_size = 3;
    smoothresults = movmean(diag(results{1}),wdw_size);
    plot(time,smoothresults,'LineWidth',1.5,'Color','k')
    yline(0.5,'r--')
    xline(0,'r')
    xlabel('Time (s)')
    ylabel('Accuracy')

    subplot(2,1,2)
    imagesc(time,time,results{1},[c_min,c_max]); axis xy; colorbar;
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
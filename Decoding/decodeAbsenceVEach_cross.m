function decodeAbsenceVEach_cross(cfg0,subject)

%% Decodes Zero vs each number separately in binary,cross-modal fashion. 
%This is to examine the discrimnability of 0 vs other numbers across domains. 
% To see if we still see a graded 0 in the cross-domain case. Should see that
% 0 is most discriminable from 5, and least from 1

%% Save Path
 
outputDir = fullfile(cfg0.root,cfg0.output_path,subject);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

%% Load MEG data
disp('loading..')
disp(subject)
dot_data = load(fullfile(cfg0.root,'EpochedData',subject,'dot_trials.mat'));
dot_data = dot_data.dot_trials;
dot_time = dot_data.time{1};
disp('loaded data')

sym_data = load(fullfile(cfg0.root,'EpochedData',subject,'sym_trials.mat'));
sym_data = sym_data.sym_trials;
sym_time = sym_data.time{1};
disp('loaded data')

det_data = load(fullfile(cfg0.root,'EpochedData',subject,'det_trials.mat'));
det_data = det_data.det_trials;
det_time = det_data.time{1};
disp('loaded data')

%% Remove No Resp Trials
cfgS = [];
cfgS.trials = dot_data.trialinfo(:,3) ~= 0 ;
dot_data = ft_selectdata(cfgS,dot_data);
cfgS.trials = sym_data.trialinfo(:,3) ~= 0 ;
sym_data = ft_selectdata(cfgS,sym_data);
cfgS.trials = det_data.trialinfo(:,3) ~= 0 ;
det_data = ft_selectdata(cfgS,det_data);


%% Get Trial x Channels x Time Matrix For Each Task
cfgS = [];
cfgS.keeptrials = true;
cfgS.channel=cfg0.channel;
dot_data = ft_timelockanalysis(cfgS,dot_data);
sym_data = ft_timelockanalysis(cfgS,sym_data);
det_data = ft_timelockanalysis(cfgS,det_data);

%% Select data
cfgS = [];
cfgS.trials = det_data.trialinfo(:,6) ~= 0.3 & det_data.trialinfo(:,6) ~= 0; % (no stim absent or stim-high contrast trials)
det_data = ft_selectdata(cfgS,det_data);

%% Balance Responses per Contrast Level
det_data = balanceContrasts(det_data);

%% Smooth Data
smoothed_dot_data = zeros(size(dot_data.trial));
for trial = 1:size(dot_data.trial,1)
    smoothed_dot_data(trial,:,:) = ft_preproc_smooth(squeeze(dot_data.trial(trial,:,:)),cfg0.nMeanS);
end
smoothed_arabic_data = zeros(size(sym_data.trial));
for trial = 1:size(sym_data.trial,1)
    smoothed_arabic_data(trial,:,:) = ft_preproc_smooth(squeeze(sym_data.trial(trial,:,:)),cfg0.nMeanS);
end
smoothed_det_data = zeros(size(det_data.trial));
for trial = 1:size(det_data.trial,1)
    smoothed_det_data(trial,:,:) = ft_preproc_smooth(squeeze(det_data.trial(trial,:,:)),cfg0.nMeanS);
end

for num = 1:5
    smoothed_arabic_data_tmp = smoothed_arabic_data(sym_data.trialinfo(:,6) == 0 | sym_data.trialinfo(:,6) == num,:,:);
    arabic_labels = sym_data.trialinfo(sym_data.trialinfo(:,6) == 0 | sym_data.trialinfo(:,6) == num,6);
    smoothed_dot_data_tmp = smoothed_dot_data(dot_data.trialinfo(:,6) == 0 | dot_data.trialinfo(:,6) == num,:,:);
    dot_labels = dot_data.trialinfo(dot_data.trialinfo(:,6) == 0 | dot_data.trialinfo(:,6) == num,6);
    

    %Fix labels
    arabic_labels = (arabic_labels == 0)+1; % 0 = 2
    dot_labels = (dot_labels == 0)+1;
    det_labels = det_data.trialinfo(:,3); %absence= 2

    cfgS = [];
    cfgS.classifier = 'lda';
    cfgS.metric = 'auc';
    cfgS.preprocess ={'undersample','average_samples'};
    cfgS.repeat = 1;
    [test_arabic_auc,~] = mv_classify_timextime(cfgS,smoothed_det_data,det_labels,smoothed_arabic_data_tmp,arabic_labels);
    [test_dot_auc,~] = mv_classify_timextime(cfgS,smoothed_det_data,det_labels,smoothed_dot_data_tmp,dot_labels);

    test_arabic_auc = diag(test_arabic_auc');%coz mvpa light has axes switched
    test_dot_auc = diag(test_dot_auc');%coz mvpa light has axes switched

    save(fullfile(outputDir,[cfg0.output_prefix{1},'0_vs_',num2str(num)]),'test_arabic_auc');
    save(fullfile(outputDir,[cfg0.output_prefix{2},'0_vs_',num2str(num)]),'test_dot_auc');

end


end

function EpochTrialsDet(cfg0,subject)


saveDir                     = fullfile(cfg0.outdir,subject);
if ~exist(saveDir,'dir'); mkdir(saveDir); end

%% Datasets 
data = load(fullfile(cfg0.datadir,subject,'det_data.mat'));
data = struct2cell(data); data = data{1};

%for subj 15
if strcmp(subject,'sub015')
    disp('Removing block 9 for sub 15')
cfgT = [];
cfgT.trials = data.trialinfo(:,1) ~=9;
data = ft_selectdata(cfgT,data);
end

%now demean 
cfgBL = []; 
cfgBL.demean = 'yes'; 
cfgBL.baselinewindow         = [-0.2 0];
data = ft_preprocessing(cfgBL,data);

%redefine trl
cfg =[];
cfg.toilim = [-cfg0.prestim cfg0.poststim];
trials = ft_redefinetrial(cfg,data);

%Remove any trials with no stim in stim zone
figure;
pdiodeIdx = find(strcmp(trials.label, 'UADC004'));
to_remove = [];
t = trials.time{1};
on = find(t == -0.06); %start check for pdiode dip here
off = find(t== 0.3); %end check here
for i = 1:length(trials.time)
    pd = trials.trial{i}(pdiodeIdx,:);
    if ~ any(abs(diff(pd(on:off))) > 0.5) %remove trial if no dip in photodiode where stim should be
        fprintf('Removing Trial %d \n',i)
        to_remove = [i to_remove];
        continue
    end
    plot(trials.time{1},trials.trial{i}(pdiodeIdx,:))
    hold on
end
xline(0,'r')
title('Det Trials')
xlabel('Time')
ylabel('Photodiode Signal')
drawnow;

%remove bad trial(s)
all_trials = 1:length(trials.trialinfo);
all_trials(to_remove) = [];
cfg = [];
cfg.trials = all_trials;
trials = ft_selectdata(cfg,trials);

det_trials  = trials; clear trials
det_trials.trialnumbers = (1:length(det_trials.trial))';
save(fullfile(saveDir,cfg0.saveName),'det_trials','-v7.3')



end
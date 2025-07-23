function EpochTrialsNum(cfg0,subject)
% epochs trials in a different way in case original script doesnt work 
% for example due to the PD signal being too noisy

saveDir                     = fullfile(cfg0.outdir,subject);
if ~exist(saveDir,'dir'); mkdir(saveDir); end

%% Datasets 
if contains(cfg0.saveName,'sym')
    data = load(fullfile(cfg0.datadir,subject,'sym_data.mat'));
    data = struct2cell(data); data = data{1};
elseif contains(cfg0.saveName,'dot')
    data = load(fullfile(cfg0.datadir,subject,'dot_data.mat'));
    data = struct2cell(data); data = data{1};
end


%% Create new trl matrix

% Get photodiode signal
cfg = [];
cfg.channel = {'UADC004'} ;
cfg.hpfilter = 'yes';
cfg.hpfreq = 1.5;
lightDiodeSignal = ft_preprocessing(cfg, data);


% Define samples prestim and poststim
trl = [];
presamples = cfg0.prestim * data.fsample;
postsamples = cfg0.poststim * data.fsample;
%figure;
for iTrial = 1:length(lightDiodeSignal.trial)
    %figure
    hold on
    lightDiodeSignal.trial{iTrial} = smoothdata(lightDiodeSignal.trial{iTrial},'gaussian',5);
    %plot(lightDiodeSignal.trial{iTrial})

    PD_on = lightDiodeSignal.trial{iTrial} < (mean(lightDiodeSignal.trial{iTrial})-(0.25*std(lightDiodeSignal.trial{iTrial}))); %get when the PD is on (when is less than 1 SD below mean)
    
    pdind = find(PD_on);
    %xline(pdind)
    [~,ind] = maxk(diff(pdind),10);
    ind = sort(pdind(ind+1))-2; %hard coded '-2' seems to align the photodiode better to 0
    %xline(ind)
    %%

    for t = 1:length(ind) %loop over the 10 stims in each trial and make their own row in the trl matrix (i.e. make them their own trial)
        stimOn = ind(t)+data.sampleinfo(iTrial,1); % stim onset in samples relative to beginning of whole experiment
        trlbegin = round(stimOn - presamples); % start 0.075 s before stim on 
        trlend   = round(stimOn + postsamples); %end 0.8s after stim on
        offset   = round(presamples); %this represents where the stim comes on within the trial samples (i.e. the offset between trial start and stim start)
        newtrl   = [trlbegin trlend -offset];
        trl      = [trl; newtrl];
        
    end
   
    
end


%redefine trl
cfg =[];
cfg.trl = trl;
trials = ft_redefinetrial(cfg,data);

%Plot new, shorter, trials and check they are all aligned with the stim at 0
pdiodeIdx = find(strcmp(trials.label, 'UADC004'));

figure;
for i = 1:length(trials.time)
    plot(trials.time{1},trials.trial{i}(pdiodeIdx,:))
    hold on
end
xline(0,'r')
if contains(cfg0.saveName,'sym')
    title('Symbolic Trials')
elseif contains(cfg0.saveName,'dot')
    title('Dot Trials')
end
xlabel('Time')
ylabel('Photodiode Signal')
drawnow;

% Baseline Correct
cfgS                        = [];
cfgS.demean                 = 'yes'; % baseline correction on 75 ms before stim 
cfgS.baselinewindow         = [-0.1 0];
num_trials          = ft_preprocessing(cfgS,trials);

%{
%Pivot the trialinfo. Columns: Run, trial number, task, numeral, colour, response, RT, correct
arabic_trials.trialinfo    = pivotTrialsArabic(num_trials.trialinfo);
% add trialnumbers for later
arabic_trials.trialnumbers = (1:length(arabic_trials.trial))';

%}

% Preallocate new matrix
newTrialInfo = zeros(size(num_trials.trialinfo, 1), 7); % First 5 columns + 2 new columns

% Copy the first 5 columns (unchanged)
newTrialInfo(:, 1:5) = num_trials.trialinfo(:, 1:5);

% Process each row to extract the relevant information
for rowIdx = 1:size(newTrialInfo, 1)
    % Determine the sequence position (1-10) based on the row index
    sequencePos = mod(rowIdx - 1, 10) + 1; % 1 to 10 repeating
    
    % Extract the number from the corresponding sequence column
    newTrialInfo(rowIdx, 6) = num_trials.trialinfo(rowIdx, 5 + sequencePos);
    
    %Colour info
    newTrialInfo(rowIdx, 7) = num_trials.trialinfo(rowIdx, 15 + sequencePos); % Store the color index
end
if contains(cfg0.saveName,'sym')
    sym_trials  = num_trials; clear num_trials trials
    sym_trials.trialinfo = newTrialInfo;
    sym_trials.trialnumbers = (1:length(sym_trials.trial))';
    save(fullfile(saveDir,cfg0.saveName),'sym_trials','-v7.3')
elseif contains(cfg0.saveName,'dot')
    dot_trials  = num_trials; clear num_trials trials
    dot_trials.trialinfo = newTrialInfo;
    dot_trials.trialnumbers = (1:length(dot_trials.trial))';
    save(fullfile(saveDir,cfg0.saveName),'dot_trials','-v7.3')
end



end
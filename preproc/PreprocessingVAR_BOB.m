function PreprocessingVAR_BOB(cfg0,subject)
% function PreprocessingVAR(cfg0,subjectdata)
% Visual artefact rejection as part of the preprocessing progress. 

% some settings
outputDir               = fullfile(cfg0.datadir,'VARData',subject);
if ~exist(outputDir,'dir'); mkdir(outputDir); end
PreprocData             = fullfile(cfg0.root,subject);

if ~exist(fullfile(outputDir,[cfg0.saveName,'.mat']),"file")
    
    
    %% Loop over data segments
    
    [~,name,~] = fileparts(cfg0.dataName);
    artSave    = ['art' name(1:3)];   
       
    % load the data
    data = load(fullfile(PreprocData,cfg0.dataName));
    data = struct2cell(data);
    data = data{1};
    
    % get some basic values
    time                    = data.time{1};
    nTrials                 = length(data.trial);
    
    % select only the MEG channels
    cfg                     = [];
    cfg.channel             = 'MEG';
    tmp_data                = ft_selectdata(cfg,data);
    tmp_data.trialinfo      = [tmp_data.trialinfo, (1:nTrials)'];
    
    % Overall artifacts
    cfg                     = [];
    cfg.method              = 'summary';
    cfg.keeptrial = 'no';
    tmp_data_overall        = ft_rejectvisual(cfg, tmp_data); % save them
    removed_n_overall       = setdiff(1:nTrials, tmp_data_overall.trialinfo(:, end));
    clear tmp_data_overall;
    
    % Muscle artifacts
    iNan = cell(nTrials); %change NaN to 0
    for in = 1:length(tmp_data.trial)
        iNan{in} = isnan(tmp_data.trial{in});
        tmp_data.trial{in}(iNan{in}) = 0;
    end
    
    cfg                     = [];
    cfg.hpfilter            = 'yes';
    cfg.hpfreq              = 100;
    tmp_data_filtered       = ft_preprocessing(cfg, tmp_data);
    
    for in = 1:length(tmp_data.trial)
        %change NaN back to NaN 
        tmp_data_filtered.trial{in}(iNan{in}) = nan;
    end
    cfg.keeptrial = 'no';
    tmp_data_muscle         = ft_rejectvisual([], tmp_data_filtered);
    removed_n_muscle        = setdiff(1:nTrials, tmp_data_muscle.trialinfo(:, end));
    clear tmp_data_muscle tmp_data_filtered;
    
    
    % Blinks during stimulus
    if cfg0.blinks
        cfg                     = [];
        cfg.channel             = {'UADC001','UADC003'}; %eye tracker X, pupil diameter
        cfg.latency             = [-0.2 0.8]; % to avoid error when trial cut off by following trial
        tmp_data                = ft_selectdata(cfg, data);
        tmp_data.trialinfo      = [tmp_data.trialinfo, (1:nTrials)'];
        disp(size(tmp_data.trial))
        disp(size(tmp_data.time{1}))

        X = reshape(cell2mat(tmp_data.trial), [2, length(tmp_data.time{1}), length(tmp_data.trial)]);
        X_EOGv = squeeze(X(1, :, :))';
        X_EOGv = (X_EOGv - mean(X_EOGv(:))) ./ std(X_EOGv(:)); %z score standardise
        X_pupDil = squeeze(X(2, :, :))';
        X_pupDil = (X_pupDil - mean(X_pupDil(:))) ./ std(X_pupDil(:)); % z score 
        
        cutoff_z_EOGv = mean(X_EOGv(:))+(2*std(X_EOGv(:))); %cut off is 2 SD +- the mean
        cutoff_z_pupDil = mean(X_pupDil(:))+(2*std(X_pupDil(:)));%cut off is 2 SD +- the mean
        
        %indexes where signal is greater than cutoff
        blinkMask               = (abs(X_EOGv) > cutoff_z_EOGv) | (abs(X_pupDil) > cutoff_z_pupDil);
        tmp                     = tmp_data.time{1} >= cfg0.stimOn(1) & tmp_data.time{1} <= cfg0.stimOn(2); %where blinks co occur with stimulus
        stimOnMask              = repmat(tmp,[nTrials,1]);
        removed_n_blinks        = tmp_data.trialinfo(any(stimOnMask & blinkMask, 2), end);
        
        clear tmp_data blinkMask stimOnMask
    end
    
    
    % Inspect potentially contaminated trials
    cfg                     = [];
    cfg.channel             = 'MEG';
    cfg.megscale = 1;
    cfg.artfctdef.overall.artifact = data.sampleinfo(removed_n_overall, :);
    cfg.artfctdef.muscle.artifact = data.sampleinfo(removed_n_muscle, :);
    if cfg0.blinks
        cfg.artfctdef.blinks.artifact = data.sampleinfo(removed_n_blinks, :);    
    end
    cfg.preproc.hpfilter    = 'no';
    cfg.preproc.hpfreq      = 100;
    cfg.renderer            = 'painters';
    cfgart                  = ft_databrowser(cfg, data);
    
    clear removed_n_overall removed_n_muscle removed_n_blinks    
    
    % Save the artifacts
    save(fullfile(outputDir,artSave),'cfgart')
    
    % Reject the artifacts
    cfgart.artfctdef.reject = 'complete';
    data                    = ft_rejectartifact(cfgart, data);
    
    % Save the data and clean up
    save(fullfile(outputDir,cfg0.saveName),'data','-v7.3')
    clear data cfgart
else
    fprintf('Already done VAR for subject %s\n ',subject)
end
        
end


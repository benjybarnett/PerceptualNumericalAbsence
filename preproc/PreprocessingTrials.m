function [cfg,data] = PreprocessingTrials(cfg0,subject)
% function PreprocessingTrials(subjectID)
% we preprocess with respect to the photodiode

%% Datasets 
raw_data_dir = fullfile(cfg0.datadir,subject,'meg','raw');

dataSets = str2fullfile(raw_data_dir,cfg0.wildcard);
disp(dataSets(1))
nDataSets = length(dataSets);
sprintf('%i data sets found',nDataSets)


%% Some settings
saveDir                     = fullfile(cfg0.root,subject);
if ~exist(saveDir,'dir'); mkdir(saveDir); end

cfgS                        = [];
cfgS.continuous             = 'yes';
cfgS.hpfilter               = 'yes';
cfgS.hpfreq                 = 0.5;
cfgS.dftfilter              = 'yes';
cfgS.demean                 = 'no'; % baseline correction on 200 ms before stim 
%cfgS.padding                = 10;

%% Numeral Task

%% Load behavioural data for trial info
%tech pilot went dots, det, symbolic
trialInfo = load(fullfile(cfg0.datadir,subject,'meg','trial_data','main','dot_data.mat'));
dotTrialInfo = trialInfo.results_dots;
nTrls = 22;
clear results_dots

% get the data per block 
dataS = cell(nDataSets,1);

for d = 1:nDataSets
    
    
    %fprintf('\t GETTING THE DATA FOR BLOCK %d OUT OF %d \n',d,nDataSets)
    cfg                         = cfgS;
    cfg.dataset                 = dataSets{d}; 

    % Do Filtering before epoching
    data                  = ft_preprocessing(cfg); 

    %cfg.trialfun                = 'trialfun_photodiadBOB';
    cfg.trialdef.eventtype      = 'UPPT001';
    cfg.trialdef.eventvalue     = cfg0.dotEventValue; % stimulus 1
    cfg.trialdef.pdiodetype     = 'UADC004';
    cfg.trialdef.prestim        = cfg0.prestimNum+0.25; %we add additional time in case PD alginment requires us to move trial's 0 point
    cfg.trialdef.poststim       = cfg0.poststimNum+0.25;
    cfg.trialdef.nTrls          = nTrls; % number of trials per block
    cfg.plot                    = cfg0.plot;
    cfg                         = ft_definetrial(cfg);
    
    % get it
    data                  = ft_redefinetrial(cfg,data);    
    
    %align with photodiode
    cfg.prestim = cfg0.prestimNum;
    cfg.poststim = cfg0.poststimNum;
    dataS{d} = AlignPDiode(cfg,data,subject,d);
    clear data
    fclose('all');
    
end


% append data
cfg = []; dot_data = ft_appenddata(cfg,dataS{:}); clear dataS

%Baseline correction
cfgBL = [];
cfgBL.demean                 = 'yes'; % baseline correction on 200 ms before stim 
cfgBL.baselinewindow         = [-cfg0.prestimNum 0];
dot_data = ft_preprocessing(cfgBL,dot_data);

% add trialnumbers for later
dot_data.trialnumbers = (1:length(dot_data.trial))';

% add trial-info
dot_data.trialinfo = dotTrialInfo;

% downsample
cfg.resamplefs              = 250;
dot_data                        = ft_resampledata(cfg, dot_data); % resample the data

% fix sample info 
dot_data = fixsampleinfo(dot_data);
dot_data = rmfield(dot_data, 'cfg');

% save and clean up
save(fullfile(saveDir,['dot_',cfg0.saveName]),'dot_data','-v7.3')
clear cfg data dot_data dotTrialInfo

%% Symbolic

%% Load behavioural data for trial info
trialInfo = load(fullfile(cfg0.datadir,subject,'meg','trial_data','main','sym_data.mat'));
symTrialInfo = trialInfo.results_symbolic;
clear results_symbolic

% get the data per block 
dataS = cell(nDataSets,1);

for d = 1:nDataSets
    
    
    %fprintf('\t GETTING THE DATA FOR BLOCK %d OUT OF %d \n',d,nDataSets)
    cfg                         = cfgS;
    cfg.dataset                 = dataSets{d}; 
    
    % Do Filtering before epoching
    data                  = ft_preprocessing(cfg); 

    %cfg.trialfun                = 'trialfun_photodiadBOB';
    cfg.trialdef.eventtype      = 'UPPT001';
    cfg.trialdef.eventvalue     = cfg0.symEventValue; % stimulus 1
    cfg.trialdef.pdiodetype     = 'UADC004';
    cfg.trialdef.prestim        = cfg0.prestimNum+0.25; %we add additional time in case PD alginment requires us to move trial's 0 point
    cfg.trialdef.poststim       = cfg0.poststimNum+0.25;
    cfg.trialdef.nTrls          = nTrls; % number of trials per block
    cfg.plot                    = cfg0.plot;
    cfg                         = ft_definetrial(cfg);
    
    % get it
    data                  = ft_redefinetrial(cfg,data);    
    
    %align with photodiode
    cfg.prestim = cfg0.prestimNum;
    cfg.poststim = cfg0.poststimNum;
    dataS{d} = AlignPDiode(cfg,data,subject,d);
    clear data
    fclose('all');
    
end


% append data
cfg = []; sym_data = ft_appenddata(cfg,dataS{:}); clear dataS

%Baseline correction
cfgBL = [];
cfgBL.demean                 = 'yes'; % baseline correction on 200 ms before stim 
cfgBL.baselinewindow         = [-cfg0.prestimNum 0];
sym_data = ft_preprocessing(cfgBL,sym_data);

% add trialnumbers for later
sym_data.trialnumbers = (1:length(sym_data.trial))';

% add trial-info
sym_data.trialinfo = symTrialInfo;

% downsample
cfg.resamplefs              = 250;
sym_data                        = ft_resampledata(cfg, sym_data); % resample the data

% fix sample info 
sym_data = fixsampleinfo(sym_data);
sym_data = rmfield(sym_data, 'cfg');

% save and clean up
save(fullfile(saveDir,['sym_',cfg0.saveName]),'sym_data','-v7.3')
clear cfg data sym_data

%% Detection

% Some settings
cfgS                        = [];
cfgS.continuous             = 'yes';
cfgS.hpfilter               = 'yes';
cfgS.hpfreq                 = 0.5;
cfgS.dftfilter              = 'yes';
cfgS.demean                 = 'no'; 
%cfgS.padding                = 10;

trialInfo = load(fullfile(cfg0.datadir,subject,'meg','trial_data','main','detection_data.mat'));
detTrialInfo = trialInfo.results_detection;
clear results_detection

% get the data per block 
dataS = cell(nDataSets,1);
nTrls = 60;
for d = 1:nDataSets
    
    
    %fprintf('\t GETTING THE DATA FOR BLOCK %d OUT OF %d \n',d,nDataSets)
    cfg                         = cfgS;
    cfg.dataset                 = dataSets{d}; 

    % Do Filtering before epoching
    data                  = ft_preprocessing(cfg); 

    % Epoching
    cfg = [];
    cfg.dataset                 = dataSets{d}; 
    cfg.trialdef.eventtype      = 'UPPT001';
    cfg.trialdef.eventvalue     = cfg0.detEventValue; % stimulus 1
    cfg.trialdef.pdiodetype     = 'UADC004';
    cfg.trialdef.prestim        = cfg0.prestimDet+0.25; %we add additional time in case PD alginment requires us to move trial's 0 point
    cfg.trialdef.poststim       = cfg0.poststimDet+0.25;
    cfg.trialdef.nTrls          = nTrls; % number of trials per block
    cfg.plot                    = cfg0.plot;
    cfg                         = ft_definetrial(cfg);
    
    % get it
    data =                  ft_redefinetrial(cfg, data);

    
    %align with photodiode
    cfg.prestim = cfg0.prestimDet;
    cfg.poststim = cfg0.poststimDet;
    dataS{d} = AlignPDiode(cfg,data,subject,d);
    clear data
    fclose('all');
    
end

% append data
cfg = []; det_data = ft_appenddata(cfg,dataS{:}); clear dataS

% add trialnumbers for later
det_data.trialnumbers = (1:length(det_data.trial))';

% add trial-info
det_data.trialinfo = detTrialInfo;

% downsample
cfg.resamplefs                  = 250;
det_data                        = ft_resampledata(cfg, det_data); % resample the data

% fix sample info 
det_data = fixsampleinfo(det_data);
det_data = rmfield(det_data, 'cfg');


% save and clean up
save(fullfile(saveDir,['det_',cfg0.saveName]),'det_data','-v7.3')
clear cfg data det_data

end

clear;
set(0,'DefaultFigureWindowStyle','docked');
%% SET WD

%% Paths
addpath('preproc')
addpath('utilities')
addpath('Decoding')
addpath('fieldtrip-master-MVPA\')
addpath('MVPA-Light-master/startup')

%% Set deaults
ft_defaults;
startup_MVPA_Light;

%% Absence Project Analysis
subjects = {
    'sub001'
    'sub002'
    'sub003'
    'sub004'
    'sub005'
    'sub006'
    'sub007'
    'sub008'
    'sub009' 
    'sub010'
    'sub011'
    'sub012'
    'sub013' 
    'sub014'
    'sub015'
    'sub016'
    'sub017' 
    'sub018' 
    'sub019'
    'sub020' 
    'sub021'
    'sub022'
    'sub023'
    'sub024'
    'sub025'
    %'sub026' %Chance level on numerical tasks
    'sub027'
    'sub028'
    'sub029'
    'sub030'
    }; 


%% Behavioural Analysis
for subj = 1:length(subjects)
    subject = subjects{subj};
    disp(subject)
    
    cfg = [];
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\Raw';
    cfg.plot = true;
    BehaviouralAnalysis(cfg,subject);

end
GroupBehaviouralAnalysis(cfg,subjects);

%% Preprocessing
for subj = 1:length(subjects)
    subject = subjects{subj};
    disp(subject)
    
    %experimental tasks
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\PreprocData';
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\Raw';
    cfg.dotEventValue = 1;
    cfg.symEventValue = 3;
    cfg.detEventValue = 5;
    cfg.prestimNum = 0.5;
    cfg.poststimNum = 4;
    cfg.prestimDet = 0.4; 
    cfg.poststimDet = 1.5; % Go quite far ahead in case we look at response-locked analyses
    cfg.saveName = 'data_preproc';
    cfg.wildcard = '*SF032*.ds';
    cfg.plot =true;
    PreprocessingTrials(cfg,subject);

end

%% VAR
for subj = 1:length(subjects)
    subject= subjects{subj};

    % Dots
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\PreprocData';
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.stimOn = [0 0]; 
    cfg.dataName = 'dot_data_preproc';
    cfg.saveName = 'dot_data_VAR';
    cfg.blinks = 0;
    PreprocessingVAR_BOB(cfg,subject);
   
    % Symbolic
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\PreprocData';
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.stimOn = [0 0.25];
    cfg.dataName = 'sym_data_preproc';
    cfg.saveName = 'sym_data_VAR';
    cfg.blinks = 0;
    PreprocessingVAR_BOB(cfg,subject);

    % Detection
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\PreprocData';
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.stimOn = [0 0.25];
    cfg.dataName = 'det_data_preproc';
    cfg.saveName = 'det_data_VAR';
    cfg.blinks = 1;
    PreprocessingVAR_BOB(cfg,subject);
    
    
end


%% Independent Components Analysis (ICA)
for subj = 1:length(subjects)
    
    subject= subjects{subj};
    disp(subject)

    tic
    % Symbolic
    cfg = [];
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.inputData = 'sym_data_VAR'; 
    cfg.compOutput = 'sym_comp';
    cfg.outputData = 'sym_data';
    PreprocessingICA_BOB(cfg,subject);
    toc

    % Dot
    cfg = [];
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.inputData = 'dot_data_VAR'; 
    cfg.compOutput = 'dot_comp';
    cfg.outputData = 'dot_data';
    PreprocessingICA_BOB(cfg,subject);
    
    % Detection
    cfg = [];
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.inputData = 'det_data_VAR'; 
    cfg.compOutput = 'det_comp';
    cfg.outputData = 'det_data';
    PreprocessingICA_BOB(cfg,subject);

end


%% Eyeball Preprocessed Data
for subj = 1:length(subjects)

    % Eyeball Dot Task
    subject = subjects{subj};
    cfg = [];
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.trials = 'all';
    cfg.channel = 'meg';
    cfg.layout = 'CTF275.lay';
    cfg.datatype = 'CleanData';
    
    cfg.file_name = 'dot_data';
    EyeballData(cfg,subject);

    % Eyeball Dymbolic Task
    cfg.file_name = 'sym_data';
    EyeballData(cfg,subject);
    
    % Eyeball Detection Task
    cfg.file_name = 'det_data';
    EyeballData(cfg,subject);
end


%% Split Into Smaller Epochs
for subj = 1:length(subjects)
    subject = subjects{subj};
    disp(subject)
    
  
    % Symbolic
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\CleanData';
    cfg.saveName = 'sym_trials.mat';
    cfg.outdir = 'D:\bbarnett\Documents\SensoryAbsence\data\EpochedData';
    cfg.poststim = 0.8;
    cfg.prestim = 0.1;
    EpochTrialsNum(cfg,subject);

    % Dots
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\CleanData';
    cfg.saveName = 'dot_trials.mat';
    cfg.outdir = 'D:\bbarnett\Documents\SensoryAbsence\data\EpochedData';
    cfg.poststim = 0.8;
    cfg.prestim = 0.1;
    EpochTrialsNum(cfg,subject);
    

    % Detection
    tic
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\CleanData';
    cfg.saveName = 'det_trials.mat';
    cfg.outdir = 'D:\bbarnett\Documents\SensoryAbsence\data\EpochedData';
    cfg.poststim = 0.8;
    cfg.prestim = 0.1;
    EpochTrialsDet(cfg,subject);
    toc

end

%% Colour Decoding
for subj = 1:length(subjects)

    subject = subjects{subj};

    % Symbolic
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryDecoding';
    cfg.output_prefix = 'symbolic_colour';
    cfg.task = 'sym';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.label = 'data.trialinfo(:,7)'; % Colour
    cfg.nAvgSamples = 10; %Number of samples to average over before decoding
    cfg.metric = 'acc';
    cfg.plot = 1;
    cfg.pltTitle = 'Colour Decoding in Symbolic Task';
    DecodeBinary(cfg,subject);

    % Dots
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryDecoding';
    cfg.output_prefix = 'dots_colour';
    cfg.task = 'dot';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.label = 'data.trialinfo(:,7)'; % Colour
    cfg.nAvgSamples = 10;
    cfg.metric = 'acc';
    cfg.plot = 1;
    cfg.pltTitle = 'Colour Decoding in Dots Task';
    DecodeBinary(cfg,subject);

    
    % Detection
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryDecoding';
    cfg.output_prefix = 'det_colour';
    cfg.task = 'det';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.label = 'data.trialinfo(:,7) '; % Colour
    cfg.nAvgSamples = 10;
    cfg.metric = 'acc';
    cfg.plot = 1;
    cfg.pltTitle = 'Colour Decoding in Detection Task';
    DecodeBinary(cfg,subject);
end


%% Detection Response Decoding
for subj = 1:length(subjects)

    subject = subjects{subj};

    % Detection Responses
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryDecoding';
    cfg.task = 'det';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.selData = 'data.trialinfo(:,6) ~= 0.3 & data.trialinfo(:,6) ~= 0'; % (no stim absent or stim-high contrast trials)
    cfg.label = 'data.trialinfo(:, 3);'; % Response: Present or Absent 
    cfg.nAvgSamples = 5;
    cfg.metric = 'acc';
    cfg.plot = 1;
    cfg.balanceContrasts = 1;
    cfg.pltTitle = 'Response Decoding in Detection Task (Balanced Contrasts)';
    cfg.output_prefix = 'det_response_CBalanced';
    DecodeBinary(cfg,subject);
end


%% Detection Stimulus Present Decoding
for subj = 1:length(subjects)

    subject = subjects{subj};

    % Detection High Contrast Hits vs. CRs
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryDecoding';
    cfg.output_prefix = 'det_presence_highC';
    cfg.task = 'det';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.selData = '(data.trialinfo(:,6) == 0.3 | data.trialinfo(:,6) == 0) & data.trialinfo(:,5) == 1'; % Select only high contrast hits or absent CRs
    cfg.label = 'data.trialinfo(:,2)'; % Stim Present or Absent 
    cfg.nAvgSamples = 5;
    cfg.metric = 'acc';
    cfg.plot = 1;
    cfg.pltTitle = 'High Contrast Hits vs. Correct Rejections';
    DecodeBinary(cfg,subject);
end



%% Number Decoding Multiclass
for subj = 1:length(subjects)

    subject = subjects{subj};

    % Symbolic
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/NumberMulticlass';
    cfg.output_prefix = 'sym_multi';
    cfg.task = 'sym';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.label = 'data.trialinfo(:,6)'; % Number
    cfg.nAvgSamples = 15;
    cfg.metric = 'acc';
    cfg.plot = 1;
    cfg.pltTitle = 'Multiclass Numbers in Symbolic';
    DecodeMulticlass(cfg,subject);

    % Dots
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/NumberMulticlass';
    cfg.output_prefix = 'dot_multi';
    cfg.task = 'dot';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.label = 'data.trialinfo(:,6)'; % Number
    cfg.nAvgSamples = 15;
    cfg.metric = 'acc';
    cfg.plot = 1;
    cfg.pltTitle = 'Multiclass Numbers in Dots';
    DecodeMulticlass(cfg,subject);
end



%% Cross-Colour
for subj = 1:length(subjects)
    subject = subjects{subj};

    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/CrossDecodingBinary';
    cfg.output_prefix = 'colour';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.label = 'data.trialinfo(:,7)'; % Colour
    cfg.nAvgSamples = 10; %Number of samples to average over before decoding
    cfg.metric = 'acc';
    cfg.plot = 1;
    cfg.pltTitle = 'Cross Decoding of Colour';
    CrossDecodeBinary3Way(cfg,subject);
end


%% Within Format Zero
for subj = 1:length(subjects)

    subject = subjects{subj};

    % Symbolic Zero
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryDecoding';
    cfg.task = 'sym';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.label = 'data.trialinfo(:, 6);'; %Number
    cfg.nAvgSamples = 5;
    cfg.metric = 'acc';
    cfg.plot = 0;
    cfg.pltTitle = 'symbolic Zero Decoding';
    cfg.output_prefix = 'zero_binary_symbolic';
    DecodeZeroBinary(cfg,subject);

    % Non-Symbolic Zero
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryDecoding';
    cfg.task = 'dot';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.label = 'data.trialinfo(:, 6);'; %Number
    cfg.nAvgSamples = 5;
    cfg.metric = 'acc';
    cfg.plot = 0;
    cfg.pltTitle = 'Non-Symbolic Zero Decoding';
    cfg.output_prefix = 'zero_binary_dot';
    DecodeZeroBinary(cfg,subject);
end



%% Graded Cross Decode of Numbers
for subj = 1:length(subjects)
    subject = subjects{subj};
    cfg  = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/CrossDecodingZero_v_Each';
    cfg.channels = 'MEG';
    cfg.nMeanS = 7;
    cfg.output_prefix =  {'train_arabic_','train_dot_'}; %must be number data first
    cfg.plot = false;
    cfg.channel = 'MEG';
    decodeZeroVEach_cross(cfg,subject);
end


%% Cross-Decoding Absence
for subj = 1:length(subjects)
    subject = subjects{subj};

    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/CrossDecodingBinary';
    cfg.output_prefix = 'absence';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.numLabel = 'data.trialinfo(:,6)'; % Number
    cfg.balanceNonZero = 1;
    cfg.selData = 'det_data.trialinfo(:,6) ~= 0.3 & det_data.trialinfo(:,6) ~= 0'; % (no stim absent or stim-high contrast trials)
    cfg.balanceContrasts = 1;
    cfg.detLabel = 'data.trialinfo(:, 3);'; % Response: Present or Absent 
    cfg.nAvgSamples = 10; %Number of samples to average over before decoding
    cfg.metric = 'acc';
    cfg.plot = 0;
    cfg.pltTitle = 'Cross Decoding of Absences';
    CrossDecodeBinary3Way(cfg,subject);
end

%% Graded Cross Decode of Absence with Numbers
for subj = 1:length(subjects)
    subject = subjects{subj};
    cfg  = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/CrossDecodingAbsence_v_Each';
    cfg.channels = 'MEG';
    cfg.nMeanS = 7;
    cfg.output_prefix =  {'test_arabic_','test_dot_'}; %must be number data first
    cfg.plot = false;
    cfg.channel = 'MEG';
    decodeAbsenceVEach_cross(cfg,subject);
end


%% Cross Decode Contrast with Empty Sets
for subj = 1:length(subjects)

    subject = subjects{subj};
    
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryCrossDecoding';
    cfg.output_prefix = 'contrast_emptySets';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.selData = '(det_data.trialinfo(:,6) == 0.3 | det_data.trialinfo(:,6) == 0)'; % Select only high contrast hits or absent CRs
    cfg.det_label = 'det_data.trialinfo(:,2)'; % Stim Present or Absent 
    cfg.num_label = 'dot_data.trialinfo(:,6)'; % Number
    cfg.nAvgSamples = 5;
    cfg.metric = 'acc';
    cfg.plot = 0;
    cfg.pltTitle = 'Grating Presence x Empty Sets';
    CrossDecodeBinary(cfg,subject);
end


%% Repreprocess and filter into different frequency bands
for subj = 1:length(subjects)
    subject = subjects{subj};
    disp(subject)
    
    %experimental tasks
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\PreprocBandData';
    cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\Raw';
    cfg.dotEventValue = 1;
    cfg.symEventValue = 3;
    cfg.detEventValue = 5;
    cfg.prestimNum = 0.5;
    cfg.poststimNum = 4;
    cfg.prestimDet = 0.4; 
    cfg.poststimDet = 1.5; % Go quite far ahead in case we look at response-locked analyses
    cfg.saveName = 'data';
    cfg.wildcard = '*SF032*.ds';
    cfg.plot =true;
    PreprocessingFilteringTrials(cfg,subject);

end


%% Split Frequency Data Into Smaller Epochs
bands = {'alpha'};
for b = 1:length(bands)
    band = bands{b};
        for subj = 1:length(subjects)
            subject = subjects{subj};
            disp(subject)
            
          
            % Symbolic
            cfg = [];
            cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
            cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\PreprocBandData';
            cfg.saveName = 'sym_trials.mat';
            cfg.outdir = 'D:\bbarnett\Documents\SensoryAbsence\data\EpochedBandData';
            cfg.poststim = 0.8;
            cfg.prestim = 0.1;
            EpochTrialsNumBand(cfg,subject,band);
        
            % Dots
            cfg = [];
            cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
            cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\PreprocBandData';
            cfg.saveName = 'dot_trials.mat';
            cfg.outdir = 'D:\bbarnett\Documents\SensoryAbsence\data\EpochedBandData';
            cfg.poststim = 0.8;
            cfg.prestim = 0.1;
            EpochTrialsNumBand(cfg,subject,band);
            
        
            % Detection
            tic
            cfg = [];
            cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
            cfg.datadir = 'D:\bbarnett\Documents\SensoryAbsence\data\PreprocBandData';
            cfg.saveName = 'det_trials.mat';
            cfg.outdir = 'D:\bbarnett\Documents\SensoryAbsence\data\EpochedBandData';
            cfg.poststim = 0.8;
            cfg.prestim = 0.1;
            EpochTrialsDetBand(cfg,subject,band);
            toc
        
        end
end

%% Zero Frequency Decoding
for subj = 1:length(subjects)

       subject = subjects{subj};
    
        % Symbolic Zero
        cfg = [];
        cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
        cfg.output_path = 'Analysis/BinaryFrequencyTimeDecoding';
        cfg.task = 'sym';
        cfg.channel = 'MEG';
        cfg.nMeanS = 7;
        cfg.label = 'data.trialinfo(:, 6);'; %Number
        cfg.nAvgSamples = 5;
        cfg.metric = 'acc';
        cfg.plot = 0;
        cfg.pltTitle = 'symbolic Zero Decoding';
        cfg.output_prefix = 'zero_binary_symbolic';
        DecodeZeroBinary_freqTime(cfg,subject,'alpha');
    
        % Non-Symbolic Zero
        cfg = [];
        cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
        cfg.output_path = 'Analysis/BinaryFrequencyTimeDecoding';
        cfg.task = 'dot';
        cfg.channel = 'MEG';
        cfg.nMeanS = 7;
        cfg.label = 'data.trialinfo(:, 6);'; %Number
        cfg.nAvgSamples = 5;
        cfg.metric = 'acc';
        cfg.plot = 0;
        cfg.pltTitle = 'Non-Symbolic Zero Decoding';
        cfg.output_prefix = 'zero_binary_dot';
        DecodeZeroBinary_freqTime(cfg,subject,'alpha');
    
end


%% Detection Frequency Decoding
for subj = 1:length(subjects)
    subject = subjects{subj};

    % Detection Responses
    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/BinaryFrequencyTimeDecoding';
    cfg.task = 'det';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.selData = 'data.trialinfo(:,6) ~= 0.3 & data.trialinfo(:,6) ~= 0'; % (no stim absent or stim-high contrast trials)
    cfg.label = 'data.trialinfo(:, 3);'; % Response: Present or Absent 
    cfg.nAvgSamples = 5;
    cfg.metric = 'acc';
    cfg.plot = 0;
    cfg.balanceContrasts = 1;
    cfg.pltTitle = 'Response Decoding in Detection Task (Balanced Contrasts)';
    cfg.output_prefix = 'det_response_CBalanced';
    FrequencyDecodeTime(cfg,subject,'alpha');
end


%% Frequency Cross Decode Absence
for subj = 1:length(subjects)
    subject = subjects{subj};

    cfg = [];
    cfg.root = 'D:\bbarnett\Documents\SensoryAbsence\data\';
    cfg.output_path = 'Analysis/FrequencyTimeCrossDecodingBinary';
    cfg.output_prefix = 'absence';
    cfg.channel = 'MEG';
    cfg.nMeanS = 7;
    cfg.numLabel = 'data.trialinfo(:,6)'; % Number
    cfg.balanceNonZero = 1;
    cfg.selData = 'det_data.trialinfo(:,6) ~= 0.3 & det_data.trialinfo(:,6) ~= 0'; % (no stim absent or stim-high contrast trials)
    cfg.balanceContrasts = 1;
    cfg.detLabel = 'data.trialinfo(:, 3);'; % Response: Present or Absent 
    cfg.nAvgSamples = 10; %Number of samples to average over before decoding
    cfg.metric = 'acc';
    cfg.plot = 0;
    cfg.pltTitle = 'Cross Decoding of Absences Frequency';
    CrossDecodeBinary3Way_freqTime(cfg,subject,'alpha');

end

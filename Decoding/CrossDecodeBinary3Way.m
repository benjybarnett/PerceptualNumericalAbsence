function CrossDecodeBinary3Way(cfg0,subject)

%% Save Path
outputDir = fullfile(cfg0.root,cfg0.output_path,subject,cfg0.output_prefix);
if ~exist(outputDir,'dir'); mkdir(outputDir); end 

%% Load MEG data
disp('loading..')
disp(subject)

data = load(fullfile(cfg0.root,'EpochedData',subject,'sym_trials.mat'));
sym_data = data.sym_trials;
time = sym_data.time{1};
disp('loaded symbolic data')
data = load(fullfile(cfg0.root,'EpochedData',subject,'dot_trials.mat'));
dot_data = data.dot_trials;
disp('loaded non-symbolic data')
data = load(fullfile(cfg0.root,'EpochedData',subject,'det_trials.mat'));
det_data = data.det_trials;
disp('loaded detection data')

%% Remove No Resp Trials
cfgS = [];
cfgS.trials = sym_data.trialinfo(:,3) ~= 0 ;
sym_data = ft_selectdata(cfgS,sym_data);
cfgS.trials = dot_data.trialinfo(:,3) ~= 0 ;
dot_data = ft_selectdata(cfgS,dot_data);
cfgS.trials = det_data.trialinfo(:,3) ~= 0 ;
det_data = ft_selectdata(cfgS,det_data);


%% Get Trial x Channels x Time Matrix 
cfgS = [];
cfgS.keeptrials = true;
cfgS.channel=cfg0.channel;
sym_data = ft_timelockanalysis(cfgS,sym_data);
dot_data = ft_timelockanalysis(cfgS,dot_data);
det_data = ft_timelockanalysis(cfgS,det_data);

%% Select data
if isfield(cfg0,'selData')
    cfgS = [];
    cfgS.trials = eval(cfg0.selData);
    det_data = ft_selectdata(cfgS,det_data);
end

%% Balance Responses per Contrast Level
if isfield(cfg0,'balanceContrasts') & cfg0.balanceContrasts
    det_data = balanceContrasts(det_data);
end

%% Smooth Data
smoothed_sym_data = zeros(size(sym_data.trial));
smoothed_dot_data = zeros(size(dot_data.trial));
smoothed_det_data = zeros(size(det_data.trial));
for trial = 1:size(sym_data.trial,1)
    smoothed_sym_data(trial,:,:) = ft_preproc_smooth(squeeze(sym_data.trial(trial,:,:)),cfg0.nMeanS);
end
for trial = 1:size(dot_data.trial,1)
    smoothed_dot_data(trial,:,:) = ft_preproc_smooth(squeeze(dot_data.trial(trial,:,:)),cfg0.nMeanS);
end
for trial = 1:size(det_data.trial,1)
    smoothed_det_data(trial,:,:) = ft_preproc_smooth(squeeze(det_data.trial(trial,:,:)),cfg0.nMeanS);
end

%% Binarise labels
sym_labels = eval(strcat('sym_',cfg0.numLabel));
sym_labels = sym_labels - min(sym_labels) + 1; %make sure labels start at 1
dot_labels = eval(strcat('dot_',cfg0.numLabel));
dot_labels = dot_labels - min(dot_labels) + 1; %make sure labels start at 1
det_labels = eval(strcat('det_',cfg0.detLabel));
det_labels = det_labels - min(det_labels) + 1; %make sure labels start at 1

%% Balance non-zero numbers
if isfield(cfg0,'balanceNonZero') & cfg0.balanceNonZero

    %Each non-zero numerosity appears the same number of times
    cfgB = [];
    cfgB.numNTClass = 5;
    [smoothed_sym_data,sym_labels] = UndersampleBinarise(cfgB,smoothed_sym_data, sym_labels, 1);
    [smoothed_dot_data,dot_labels] = UndersampleBinarise(cfgB,smoothed_dot_data, dot_labels, 1);

end

%% Turn all numbers to zero and not-zero
sym_labels = (sym_labels == 1)+1; %zero = 2
dot_labels = (dot_labels == 1)+1; %empty-set = 2
%det_labels = (det_labels == 1)+1; %absence = 2
%detection already has 2 = absence

%% Within Time x Time Decoding
cfgS = [];
cfgS.classifier = 'lda';
cfgS.metric = cfg0.metric;
cfgS.preprocess ={'undersample','average_samples'};
cfgS.preprocess_param{2}.group_size = cfg0.nAvgSamples;
cfgS.repeat = 1;

dataSets = {smoothed_sym_data; smoothed_dot_data; smoothed_det_data};
dataLabels = {sym_labels; dot_labels; det_labels};
taskNames = {'sym';'dot';'det'};
results = cell(length(dataSets),length(dataSets));
for train = 1:length(dataSets)
    for test = 1:length(dataSets)
        if train ~= test & (train == 3 || test == 3)
            [output,~] = mv_classify_timextime(cfgS,dataSets{train},dataLabels{train},dataSets{test},dataLabels{test});
            %Save
            results{train,test} = output;
            save(fullfile(outputDir, [taskNames{train} '_' taskNames{test} '_' cfg0.metric]),'output');
        end
    end
end


%% Plot
if cfg0.plot

    % Define row and column labels:
    rowLabels = {'Symbolic','Non-Symbolic','Detection'};
    colLabels = {'Symbolic','Non-Symbolic','Detection'};

    % 1) Compute the global min and max across all off-diagonal matrices
    allVals = [];
    for train = 1:3
        for test = 1:3
            if train ~= test && ~isempty(results{train,test})
                allVals = [allVals; results{train,test}(:)];
            end
        end
    end
    c_min = min(allVals) - 0.05;
    c_max = max(allVals) + 0.05;

    % 2) Create a figure and set up the overall title
    figure('Units','normalized','OuterPosition',[0 0 1 1]);
    sgtitle([subject, ' ', cfg0.pltTitle])  % Keep or change as needed

    % 3) Loop over the 3x3 cells and plot off-diagonal entries
    for train = 1:3
        for test = 1:3
            % Skip diagonal
            if train ~= test && ~isempty(results{train,test})
                % Pick subplot in row train, column test
                subplot(3,3,(train-1)*3 + test);
                disp((train-1)*3 + test)

                % Extract the matrix
                thisData = results{train,test};

                % Plot with global color limits
                imagesc(time, time, thisData, [c_min, c_max]);
                axis xy; axis square;
                hold on;

                % (Optional) Mark lines for time=0 and time=1
                plot([0 0], ylim, 'r');
                plot(xlim, [0 0], 'r');
                plot([1 1], ylim, 'k');
                plot(xlim, [1 1], 'k');

                % Standard axis labels for reference
                xlabel('Train Time (s)');
                ylabel('Test Time (s)');

                % Label the row (y-axis) and column (top title)
                % Only put the row label at the leftmost column of each row:
                if test == 1
                    ylabel({rowLabels{train}, 'Test Time (s)'});
                    % Combining row label + existing 'Test Time'
                    % or just use rowLabels{i}.
                end
                % Only put the column label at the top row:
                if train == 1
                    title(colLabels{test});
                end
            end
        end
    end

    % 4) Single colormap and colorbar for all subplots
    cmap = flipud(ft_colormap('RdBu'));
    colormap(cmap);

    % Colorbar (adjust position as needed)
    hcb = colorbar('Position',[0.93 0.11 0.02 0.815]);
    caxis([c_min c_max]);

    drawnow;



end


end
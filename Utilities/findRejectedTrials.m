
for subj = 1:length(subjects)
    subject = subjects{subj};


    root = 'D:\bbarnett\Documents\SensoryAbsence\data\';

    %% Detection
    art_reject = load(fullfile(root,'VARdata',subject,'artdet'));
    art_reject = struct2cell(art_reject); art_reject = art_reject{1};

    nonVAR = load(fullfile(root,'PreprocData',subject,'det_data_preproc.mat'));
    nonVAR = struct2cell(nonVAR); nonVAR = nonVAR{1};
    all_samples = nonVAR.sampleinfo;


    all_trials = art_reject.artfctdef;
    all_trials = struct2cell(all_trials);

    rej_trials = [];

    for type = 1:length(all_trials)
        trials = all_trials{type}.artifact;
        [~,rej] = ismember(trials,all_samples); 
        rej_trials = [rej(:,1); rej_trials];
        disp(rej)
        disp(trials(:,2)-trials(:,1))

    end
    rej_trials = sort(rej_trials);

    save_to = fullfile(root,'rejectedTrials',subject);
    mkdir(save_to);cd(save_to);
    save('rej_trials_det','rej_trials');

     %% Symbolic
    art_reject = load(fullfile(root,'VARdata',subject,'artsym'));
    art_reject = struct2cell(art_reject); art_reject = art_reject{1};

    nonVAR = load(fullfile(root,'PreprocData',subject,'sym_data_preproc.mat'));
    nonVAR = struct2cell(nonVAR); nonVAR = nonVAR{1};
    all_samples = nonVAR.sampleinfo;


    all_trials = art_reject.artfctdef;
    all_trials = struct2cell(all_trials);

    rej_trials = [];

    for type = 1:length(all_trials)
        trials = all_trials{type}.artifact;
        [~,rej] = ismember(trials,all_samples); 
        rej_trials = [rej(:,1); rej_trials];


    end
    rej_trials = sort(rej_trials);

    save_to = fullfile(root,'rejectedTrials',subject);
    mkdir(save_to);cd(save_to);
    save('rej_trials_sym','rej_trials');

     %% Dots
    art_reject = load(fullfile(root,'VARdata',subject,'artdot'));
    art_reject = struct2cell(art_reject); art_reject = art_reject{1};

    nonVAR = load(fullfile(root,'PreprocData',subject,'dot_data_preproc.mat'));
    nonVAR = struct2cell(nonVAR); nonVAR = nonVAR{1};
    all_samples = nonVAR.sampleinfo;


    all_trials = art_reject.artfctdef;
    all_trials = struct2cell(all_trials);

    rej_trials = [];

    for type = 1:length(all_trials)
        trials = all_trials{type}.artifact;
        [~,rej] = ismember(trials,all_samples); 
        rej_trials = [rej(:,1); rej_trials];


    end
    rej_trials = sort(rej_trials);

    save_to = fullfile(root,'rejectedTrials',subject);
    mkdir(save_to);cd(save_to);
    save('rej_trials_dot','rej_trials');

end
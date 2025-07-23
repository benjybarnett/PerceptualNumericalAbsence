function plot_BF(cfg,subjects)
% Get subject directories
subject_dirs = dir(fullfile(cfg.root,cfg.dir, '*'));
subject_dirs = subject_dirs([subject_dirs.isdir]);
subject_dirs = subject_dirs(~ismember({subject_dirs.name}, {'.', '..'}));

% Initialize variables
all_results = [];
time = load('time.mat');time = time.time';
missing_subjs = [];
for i = 1:length(subject_dirs)
    subj_dir = fullfile(cfg.root, cfg.dir, subject_dirs(i).name);
    if any(strcmp(subjects,subject_dirs(i).name))
        data_file = fullfile(subj_dir, cfg.filename);
        %disp(subj_dir)
        if exist(data_file, 'file')
            data= load(data_file); data = struct2cell(data); data = data{1};
            all_results(:,:,i) = data;
            all_diags(:,i) = diag(data);
        else
            warning('Missing results.mat in %s', subj_dir);
            missing_subjs = [missing_subjs i];
        end
    else
        warning('Skipping Bad Subject %s', subject_dirs(i).name)
        missing_subjs = [missing_subjs i];
    end
end

%Remove missing subjects from results
all_results(:,:,missing_subjs) = [];
all_diags(:,missing_subjs) = [];
disp(size(all_diags))
disp(size(all_results))


tic
[bf, bf_comp] = bayesfactor_diagonal(all_diags,'nullvalue', 0.5,'interval', [0.5 Inf]);
toc


subplot(1,2,2); plot(time,log10(abs(bf)));
yline(1,'r--')
yline(-1,'r--')



end

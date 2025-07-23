function plot_avg_decoding(cfg,subjects)
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

% Compute average results across subjects
avg_results = mean(all_results, 3);
mean_diag = mean(all_diags,2);

% Compute 95% confidence interval
num_subjects = size(all_diags, 2);
std_diag = std(squeeze(all_diags(:, :, :)), 0, 2);
ci_diag = 1.96 * (std_diag ./ sqrt(num_subjects));

% Determine color scale limits
c_min = cfg.clim(1);
c_max = cfg.clim(2);

% Plot the results
figure('units', 'normalized', 'outerposition', [0 0 1 1]);

% Smoothed diagonal plot with confidence interval
subplot(2,1,1);
wdw_size = 3;
smooth_avg_results = movmean(mean_diag, wdw_size);
smooth_ci_upper = movmean(mean_diag + ci_diag, wdw_size);
smooth_ci_lower = movmean(mean_diag - ci_diag, wdw_size);

fill([time; flipud(time)], [smooth_ci_upper; flipud(smooth_ci_lower)], [1 0.5 0.2],'EdgeColor','none','FaceAlpha',0.6);
hold on;
plot(time, smooth_avg_results, 'LineWidth', 1.5, 'Color', 'k');
yline(cfg.chance, 'r--');
xlim([min(time),max(time)]);
ylim(cfg.ylim);
xline(0, 'r');
xlabel('Time (s)', 'FontName', 'Arial', 'FontSize', 16);
ylabel('Accuracy', 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'FontName', 'Arial', 'FontSize', 14);

% Heatmap plot
subplot(2,1,2);
imagesc(time, time, avg_results, [c_min, c_max]); axis xy; cb = colorbar; cbpos = get(cb, 'Position');
xlabel('Train Time (s)', 'FontName', 'Arial', 'FontSize', 16);
ylabel('Test Time (s)', 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'FontName', 'Arial', 'FontSize', 14);
ylabel(cb, 'Accuracy', 'Rotation', 270);
set(get(cb, 'YLabel'), 'Position', [cbpos(1) + 3.25, cbpos(2) + cbpos(4)+0.1, 0],'FontName', 'Arial', 'FontSize', 16);
hold on; plot([0 0], ylim, 'r'); hold on; plot(xlim, [0 0], 'r');
hold on; plot([1 1], ylim, 'k'); hold on; plot(xlim, [1 1], 'k');
cmap = flipud(ft_colormap('RdBu'));
colormap(cmap);
axis square;
sgtitle(cfg.pltTitle, 'FontName', 'Arial', 'FontSize', 16);
drawnow;

if cfg.sigTest
    cfgP = [];
    cfgP.indiv_pval = 0.05;
    cfgP.cluster_pval = 0.05;
    cfgP.iterations = 1000;
    cfgP.paired = false;
    cfgP.tail= 'one';
    clusterPvals = cluster_based_permutationND(permute(all_results,[3,1,2]),cfg.chance,cfgP);
    disp(min(clusterPvals,[],'all'))
    sigMask = clusterPvals < 0.05;
    if any(sigMask(:))
        disp('Significant Clusters!')
    end
    % Find significant clusters
    hold on;
    significant_clusters = bwboundaries(sigMask); % Extract boundaries
    for k = 1:length(significant_clusters)
        boundary = significant_clusters{k};
        if length(boundary) > cfg.plotClustLim
            plot(time(boundary(:,2)), time(boundary(:,1)), 'k', 'LineWidth', 1.5); % Overlay cluster boundaries
        end
    end
end

end

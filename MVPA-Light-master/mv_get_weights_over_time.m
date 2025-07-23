function cf = mv_get_weights_over_time(cfg, X, clabel)
% Classification across time. A classifier is trained and validate for
% different time points in the dataset X. Cross-validation should be used
% to get a realistic estimate of classification performance.
%
% Usage:
% [perf, res] = mv_classify_across_time(cfg,X,clabel)
%
%Parameters:
% X              - [samples x features x time points] data matrix -OR-
%                  [samples x samples  x time points] kernel matrices
% clabel         - [samples x 1] vector of class labels
%
% cfg          - struct with hyperparameters:
% .classifier   - name of classifier, needs to have according train_ and test_
%                 functions (default 'lda')
% .hyperparameter - struct with parameters passed on to the classifier train
%                 function (default [])
% .metric       - classifier performance metric, default 'accuracy'. See
%                 mv_classifier_performance. If set to [] or 'none', the 
%                 raw classifier output (labels, dvals or probabilities 
%                 depending on cfg.output_type) for each sample is returned. 
%                 Use cell array to specify multiple metrics (eg
%                 {'accuracy' 'auc'}
% .time         - indices of time points (by default all time
%                 points in X are used)
% .feedback     - print feedback on the console (default 1)
%
% CROSS-VALIDATION parameters:
% .cv           - perform cross-validation, can be set to 'kfold',
%                 'leaveout', 'holdout', 'predefined' or 'none' (default 'kfold')
% .k            - number of folds in k-fold cross-validation (default 5)
% .p            - if cv is 'holdout', p is the fraction of test samples
%                 (default 0.1)
% .stratify     - if 1, the class proportions are approximately preserved
%                 in each fold (default 1)
% .repeat       - number of times the cross-validation is repeated with new
%                 randomly assigned folds (default 1)
% .fold         - if cv='predefined', fold is a vector of length
%                 #samples that specifies the fold each sample belongs to
%
% PREPROCESSING parameters:
% .preprocess         - cell array containing the preprocessing pipeline. The
%                       pipeline is applied in chronological order (default {})
% .preprocess_param   - cell array of preprocessing parameter structs for each
%                       function. Length of preprocess_param must match length
%                       of preprocess (default {})
%
% Returns:
% perf          - [time x 1] vector of classifier performances. If
%                 metric='none', perf is a [r x k x t] cell array of
%                 classifier outputs, where each cell corresponds to a test
%                 set, k is the number of folds, r is the number of 
%                 repetitions, and t is the number of time points. If
%                 multiple metrics are requested, perf is a cell array
% result        - struct with fields describing the classification result.
%                 Can be used as input to mv_statistics and mv_plot_result
% testlabel     - [r x k] cell array of test labels. Can be useful if
%                 metric='none'
%
% Note: For time x time generalisation, use mv_classify_timextime

% (c) Matthias Treder

X = double(X);
if ndims(X)~= 3, error('X must be 3-dimensional'), end

mv_set_default(cfg,'classifier','lda');
mv_set_default(cfg,'hyperparameter',[]);
mv_set_default(cfg,'metric','accuracy');
mv_set_default(cfg,'time',1:size(X,3));
mv_set_default(cfg,'feedback',1);

mv_set_default(cfg,'preprocess',{});
mv_set_default(cfg,'preprocess_param',{});

[cfg, clabel, n_classes, n_metrics] = mv_check_inputs(cfg, X, clabel);

ntime = numel(cfg.time);

% Number of samples in the classes
n = arrayfun( @(c) sum(clabel==c) , 1:n_classes);

% indicates whether the data represents kernel matrices
mv_set_default(cfg,'is_kernel_matrix', isfield(cfg.hyperparameter,'kernel') && strcmp(cfg.hyperparameter.kernel,'precomputed'));
if cfg.is_kernel_matrix,  mv_set_default(cfg,'dimension_names',{'samples','samples','time points'});
else,                     mv_set_default(cfg,'dimension_names',{'samples','features','time points'}); end

%% Get train and test functions
train_fun = eval(['@train_' cfg.classifier]);
test_fun = eval(['@test_' cfg.classifier]);

%% Train decoder over time to get weights per time point
if cfg.feedback, mv_print_classification_info(cfg,X,clabel); end

% No cross-validation, just train and test once for each
% training/testing time. This gives the classification performance for
% the training set, but it may lead to overfitting and thus to an
% artifically inflated performance.

% Initialize classifier outputs
cf_output = cell(1, 1, ntime);

% Preprocess train/test data
[~, X, clabel] = mv_preprocess(cfg, X, clabel);

for tt=1:ntime          % ---- Train and test time ----
    % Train and test data
    Xtraintest= squeeze(X(:,:,cfg.time(tt)));

    % Train classifier
    cf{tt}= train_fun(cfg.hyperparameter, Xtraintest, clabel);
    
    % Obtain classifier output (class labels or dvals)
    %cf_output{1,1,tt} = mv_get_classifier_output(cfg.output_type, cf, test_fun, Xtraintest);
end


end
function new_data = BalanceNonZeroClasses(data, labels)
% BalanceNonZeroClasses
% Balances trial data so that:
% - Trials with label 1 are equal in number to combined balanced non-1 trials
% - Each non-1 category (2–6) has an equal number of trials
%
% Inputs:
%   data   - FieldTrip data structure
%   labels - Vector of labels (e.g., data.trialinfo(:,6))
%
% Output:
%   new_data - Balanced FieldTrip data structure

    % Find trial indices for each condition
    idx1 = find(labels == 1);
    numIdx1 = numel(idx1);

    % Find indices for each of the not-1 categories (2 to 6)
    not1_categories = 2:6;
    idxEachNot1 = arrayfun(@(x) find(labels == x), not1_categories, 'UniformOutput', false);

    % Max number of trials per non-1 category based on available '1' trials
    maxPerCategory = floor(numIdx1 / 5);  % since we have 5 categories (2–6)

    % Also respect available trials in each non-1 category
    trueMinPerCategory = min(cellfun(@numel, idxEachNot1));
    countPerCategory = min(maxPerCategory, trueMinPerCategory);

    % Randomly sample from each non-1 category
    balancedIdxNot1 = cellfun(@(x) randsample(x, countPerCategory), idxEachNot1, 'UniformOutput', false);
    balancedIdxNot1 = vertcat(balancedIdxNot1{:});

    % Randomly select matching number of '1' trials
    balancedIdx1 = randsample(idx1, countPerCategory * 5);

    % Combine and sort
    selectedIdx = sort([balancedIdx1; balancedIdxNot1]);

    % Select balanced trials
    cfg = [];
    cfg.trials = selectedIdx;
    new_data = ft_selectdata(cfg, data);
end

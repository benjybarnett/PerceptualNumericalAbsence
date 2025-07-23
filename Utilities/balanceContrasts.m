function dataBalanced = balanceContrasts(data)
% BALANCECONTRASTS Balances the number of trials marked as '1' and '0' 
% in column 3 of data.trialMatrix for each unique contrast in column 6,
% by randomly selecting the smaller group size from each condition.
%
% INPUT:
%   data          - FieldTrip-like data struct with data.trialinfo
%                   (N x M matrix, where:
%                     - column 3: contains the participant response (1: present, 2: absent)
%                     - column 6: contains contrast levels)
%
% OUTPUT:
%   dataBalanced  - same data struct as input but with trials balanced 
%                   and selected using ft_selectdata
%
% EXAMPLE USAGE:
%   dataBalanced = balanceContrasts(data);
%
% NOTE: This function assumes you have FieldTrip's ft_selectdata in your path.

    % Extract the trial matrix from the data struct
    trialMatrix = data.trialinfo;

    % Identify unique contrast values in the 6th column
    contrasts = unique(trialMatrix(:,6));

    % Initialize a logical index of which trials to keep
    keepIdx = false(size(trialMatrix,1),1);

    % Loop over each contrast level
    for iC = 1:length(contrasts)
        cVal = contrasts(iC);

        % Find trials with this contrast value
        theseTrials = find(trialMatrix(:,6) == cVal);

        % Separate these by column-3 == present response (1) or absent response (2)
        presIdx = theseTrials(trialMatrix(theseTrials,3) == 1);
        absIdx  = theseTrials(trialMatrix(theseTrials,3) == 2);

        % Determine how many we can keep (match the smaller group)
        nPres = length(presIdx);
        nAbs  = length(absIdx);
        nMin  = min(nPres, nAbs);

        % Randomly permute and select nMin trials from each group
        presRandIdx = presIdx(randperm(nPres, nMin));
        absRandIdx  = absIdx(randperm(nAbs, nMin));

        % Mark these selected trials as keepers
        keepIdx(presRandIdx) = true;
        keepIdx(absRandIdx)  = true;

        fprintf('For %d-th contrast level (contrast=%g), selected %d trials from each condition.\n', ...
                iC, cVal, nMin);
    end

    fprintf('Keeping %d trials in total.\n', sum(keepIdx));

    % Use ft_selectdata to subset the data
    cfg = [];
    cfg.trials = keepIdx; 
    dataBalanced = ft_selectdata(cfg, data);

end

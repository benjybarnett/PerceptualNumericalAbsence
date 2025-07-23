function [X,Y] = getXYfromdataStruct(data,triggerVals)
% function [X,Y] = getXYfromdataStruct(data,triggerVals)
%
% divides the data into a trials by channels by samples matrix and gives
% the corresponding labels based on the triggerVals, which should be a
% nCond by 1 cell array of trigger values

% get the dimension sizes
nChannels        = size(data.trial{1,1},1);
nTrials          = numel(data.trial);
nSamples         = size(data.trial{1,1},2);

% get the X in the right format
X     = permute(reshape(cat(1,data.trial{:}),[nChannels,nTrials,nSamples]),[2,1,3]); 

% get the labels
Y = zeros(nTrials,1);
for c = 1:numel(triggerVals)    
    indices = ismember(data.trialinfo,triggerVals{c});
    Y(indices) = c;
end
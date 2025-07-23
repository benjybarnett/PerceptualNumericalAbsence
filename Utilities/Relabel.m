function [labels] = Relabel(labels)
%relabel labels so they go from 1:N
u = unique(labels);
n_classes = length(u);
if ~all(ismember(labels,1:n_classes))
    warning('Class labels should consist of integers 1 (class 1), 2 (class 2), 3 (class 3) and so on. Relabelling them accordingly.');
    newlabel = nan(numel(labels), 1);
    for i = 1:n_classes
        newlabel(labels==u(i)) = i; % set to 1:nth classes
    end
    labels = newlabel;
end
end